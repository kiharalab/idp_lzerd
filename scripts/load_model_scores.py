#!/usr/bin/env python

# Copyright (C) 2016-2017 Lenna X. Peterson, Daisuke Kihara, and Purdue University
# This file is part of IDP-LZerD.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division

import argparse
import functools
import inspect
import json
import logging
import os

from Bio import PDB
import numpy as np
import pandas as pd
pd.set_option("max_colwidth", 100)
from scipy import spatial

script_dir = os.path.abspath(os.path.dirname(inspect.getfile(inspect.currentframe())))

import shared


model_schema = """
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS allmodel
    (
        modelid INTEGER PRIMARY KEY NOT NULL,
        modelindex INTEGER NOT NULL,
        dfire REAL,
        itscore REAL NOT NULL,
        windowindex INTEGER NOT NULL,
        fragmentindex INTEGER NOT NULL,
        timestamp DATE DEFAULT (datetime('now', 'localtime')),
        FOREIGN KEY(windowindex) REFERENCES fragment(windowindex),
        FOREIGN KEY(fragmentindex) REFERENCES fragment(fragmentindex)
    );
"""


class LoadModelScoresError(RuntimeError):
    """Exception for class LoadModelScores"""
    def __init__(self, message):
        super(LoadModelScoresError, self).__init__(message)


class LoadModelScores(object):

    def __init__(self, working_dir, **kwargs):
        """CONSTRUCTOR"""

        self.modeldist_sql_data['db_name_fmt'] = os.path.join(working_dir, "{pdbid}_modeldist{windowid}.db")

        self.load(**kwargs)
        self.choose(**kwargs)

    def load(self, complexdb, complexname, receptor_chain, ligand_chain, **kwargs):
        """
        Load model scores.
        """

        fragmentselect = """
        SELECT windowindex, fragmentindex, window_wd
        FROM fragment
        JOIN window USING(windowindex)"""

        model_ct_select = "SELECT COUNT(*) FROM allmodel WHERE fragmentindex=:fragmentindex AND windowindex=:windowindex"
        modelcolumns = ["modelindex", "dfire", "itscore", "fragmentindex", "windowindex"]
        modelinsert = shared.create_insert_statement("allmodel", modelcolumns)

        pdbid = "{0}{1}{2}".format(complexname, receptor_chain, ligand_chain)
        logging.debug(pdbid)

        receptor_chain = receptor_chain.lower()[:1]
        ligand_chain = ligand_chain.lower()

        with shared.write_conn(complexdb) as conn:
            curs = conn.cursor()
            curs.execute(model_schema)
            ## INSERT MODELS
            fragment_rows = shared.conn_to_pandas(fragmentselect, conn)
            fragment_rows = fragment_rows.to_dict("records")
            for frag_row in fragment_rows:
                windowindex = frag_row['windowindex']
                fragmentindex = frag_row['fragmentindex']
                frag_wd = os.path.join(frag_row['window_wd'], str(fragmentindex))

                itscorefile = os.path.join(frag_wd, shared.itscore_file)
                model_itscores = shared.read_itscore(itscorefile, kind="model")

                nmodels = len(model_itscores)
                nmodelrows = curs.execute(model_ct_select, dict(fragmentindex=fragmentindex, windowindex=windowindex)).fetchone()[0]
                if nmodelrows < nmodels:
                    # Drop partial table
                    if nmodelrows:
                        logging.debug("Dropping partial model table")
                        curs.execute("DELETE FROM model")

                    nitscore = len(model_itscores)

                    dfirefile = os.path.join(frag_wd, shared.goap_file)
                    model_dfires = shared.read_dfire_from_goap(dfirefile, kind="model")

                    modelrows = pd.merge(model_itscores, model_dfires, on="modelindex", how="left")
                    ndi = len(modelrows)

                    if len(set((nitscore, ndi))) != 1:
                        logging.error("ITScores: %s", nitscore)
                        logging.error("IT lj dfire: %s", ndi)
                        raise LoadModelScoresError("Score number mismatch")

                    modelrows['fragmentindex'] = fragmentindex
                    modelrows['windowindex'] = windowindex
                    curs.executemany(modelinsert, modelrows.to_dict("records"))

    def choose(self, complexdb, complexname, receptor_chain, ligand_chain, **kwargs):
        """
        Choose models and load coordinates.
        """
        models_per_window = 4500

        parser = PDB.PDBParser(QUIET=True)

        scores = [dict(columns="itscore", ascending=True),
                  dict(columns="dfire", ascending=True)]

        sql_dict = self.prepare_sql(complexdb)
        n_model_q_fmt = sql_dict['n_model_q_fmt']
        window_data_q_fmt = sql_dict['window_data_q_fmt']
        model_insert = sql_dict['model_insert']

        window_q = "SELECT windowindex, window_wd, res_end - res_start + 1 AS length FROM window ORDER BY res_start"
        window_df = shared.db_to_pandas(window_q, complexdb)
        window_rows = window_df.to_dict("records")

        for x, window_row in enumerate(window_rows):
            window_wd = window_row['window_wd']
            length = window_row['length']
            n_chosen_q = n_model_q_fmt.format(**window_row)
            with shared.ro_conn(complexdb) as conn:
                n_chosen = list(conn.cursor().execute(n_chosen_q))[0][0]
            if n_chosen == models_per_window:
                continue
            get_coords = functools.partial(self.get_coords, parser=parser, chain=ligand_chain, length=length)
            window_data_q = window_data_q_fmt.format(**window_row)
            windowrows = shared.db_to_pandas(window_data_q, complexdb)
            if windowrows.empty:
                raise LoadModelScoresError("No rows for windowindex {windowindex}".format(**window_row))
            # Scale scores and create di column
            windowrows = self.scale_scores(windowrows, scores)
            # Sort by di
            windowrows = windowrows.sort_values("di")
            # Get top N rows by di
            windowrows = windowrows.head(models_per_window) 
            # Add path
            path_fmt = os.path.join(window_wd, "{fragmentindex:.0f}", "decoys", "model{modelindex:.0f}.pdb")
            windowrows['path'] = windowrows.apply(lambda x: path_fmt.format(**x), 1)
            # Prepare for insertion
            insert_rows = windowrows.apply(get_coords, 1)
            with shared.write_conn(complexdb) as conn:
                conn.cursor().executemany(model_insert, insert_rows.to_dict("records"))

        pdbid = "{0}{1}{2}".format(complexname, receptor_chain, ligand_chain)
        self.load_modeldist(complexdb=complexdb, pdbid=pdbid, window_rows=window_rows, **self.params)

    def load_modeldist(self, complexdb, pdbid, window_rows,
                       backbone_atoms, cb, stored_atoms, n_atoms, ca_index,
                       atom_overlap, res_overlap, bb_threshold, max_clash,
                       sticky_max_all, sticky_max_any, min_cosine,
                       min_ifourdist, max_ifourdist,
                       min_isixdist, max_isixdist,
                       min_itwelvedist, **kwargs):

        model_q = """
        SELECT windowindex, modelid, coordinates
        FROM model
        JOIN allmodel USING(modelid)
        JOIN fragment USING(windowindex, fragmentindex)
        ORDER BY windowindex, modelindex
        """
        models = shared.db_to_pandas(model_q, complexdb).to_dict("records")

        if not models:
            raise LoadModelScoresError("No models")

        modeldistcoords = {row['modelid']: json.loads(row['coordinates'])
                           for row in models}

        # Store consecutive pairs of window indices
        window_indices = [row['windowindex'] for row in window_rows]

        def window_dist(w1, w2):
            return abs(window_indices.index(w2) - window_indices.index(w1))

        w_models = dict()
        for w in window_indices:
            w_models[w] = [row for row in models if int(row['windowindex']) == w]

        for key, val in w_models.iteritems():
            print key, len(val)

        # Step forward by 4 residues and forward to CA atom
        mp_idx = n_atoms * (res_overlap + 1) + ca_index
        # Step backwards by 4 residues and forward to CA atom
        i_idx = -n_atoms * (res_overlap + 1) + ca_index
        # Step forward by 3 residues and forward to CA atom
        i4_idx = n_atoms * res_overlap + ca_index

        def get_mp_dist(separation):
            if separation < 1:
                raise LoadModelScoresError("Window separation %s", separation)
            elif separation == 1:
                min_dist = min_isixdist
            else:
                min_dist = min_itwelvedist
            max_dist = separation * max_isixdist
            return min_dist, max_dist

        # Pre-create np arrays of midpoint coordinates
        mp_coord_dict = dict()
        for x, window in enumerate(window_indices):
            window_length = window_rows[x]['length']
            if window_length == 4:
                # CA of 4th residue
                window_mp_idx = -4
            elif window_length > 4:
                window_mp_idx = mp_idx
            else:
                raise LoadModelScoresError("Window {windowindex} has invalid length {length}".format(**window_rows[x]))
            mp_coord_dict[window] = np.array([modeldistcoords[id['modelid']][window_mp_idx]
                                             for id in w_models[window]])

        # Calculate modeldists
        def calculate_modeldist(w1, w2):
            logging.debug("Starting %s %s", w1, w2)
            separation = window_dist(w1, w2)
            neighbors = separation == 1
            mp_min, mp_max = get_mp_dist(separation)

            # Calculate all midpoint distances simultaneously
            w1_array = mp_coord_dict[w1]
            w2_array = mp_coord_dict[w2]
            pairwise_mp_dist = spatial.distance.cdist(w1_array, w2_array)
            allowed = (pairwise_mp_dist > mp_min) & (pairwise_mp_dist < mp_max)
            it = np.nditer(allowed, flags=['multi_index'])
            while not it.finished:
                # Convert from array indices to local lists
                mp_index = it.multi_index
                mp_dist = pairwise_mp_dist[mp_index]
                i, j = mp_index
                w1model = w_models[w1][i]
                w2model = w_models[w2][j]

                id1 = w1model['modelid']
                id2 = w2model['modelid']

                row = dict(modela=id1,
                           modelb=id2,
                           mpdist=mp_dist)
                # If mp_dist constraint is met, check overlap etc.
                if it[0]:
                    skip = False
                    w1_coords = modeldistcoords[id1]
                    w2_coords = modeldistcoords[id2]

                    if neighbors:
                        # CHECK i to i+4 CA DISTANCE
                        i_ca = w1_coords[i_idx]
                        i4_ca = w2_coords[i4_idx]
                        # spatial can handle lists
                        ifourdist = spatial.distance.euclidean(i_ca, i4_ca)

                        if ifourdist < min_ifourdist or ifourdist > max_ifourdist:
                            skip = True
                        else:
                            row['ifourdist'] = ifourdist

                    # CHECK FOR BACKBONE CLASH
                    if not skip and mp_dist <= max_isixdist:
                        w2_nooverlap = w2_coords
                        if neighbors:
                            w2_nooverlap = w2_nooverlap[atom_overlap:]
                        pairwise_overlap = spatial.distance.cdist(w1_coords, w2_nooverlap)
                        clash = pairwise_overlap <= bb_threshold
                        # Number of atoms in clash
                        nclash = np.count_nonzero(clash)

                        if nclash > max_clash:
                            skip = True
                    else:
                        nclash = 0

                    if not skip and neighbors:
                        # PRECOMPUTE EDGESCORE
                        w1_sticky = w1_coords[-atom_overlap:]
                        w2_sticky = w2_coords[:atom_overlap]
                        distances = np.diag(spatial.distance.cdist(w1_sticky, w2_sticky))
                        if any(distances > sticky_max_any):
                            skip = True
                        elif all(distances > sticky_max_all):
                            skip = True
                        else:
                            row['edgescore'] = np.mean(np.square(distances))

                    if not skip and neighbors:
                        # CALCULATE COSINE

                        # N of third to last residue
                        w1_start = np.array(w1_coords[-12])
                        # C of last residue
                        w1_end = np.array(w1_coords[-2])
                        w1_vec = w1_end - w1_start
                        w1_mag = np.linalg.norm(w1_vec)

                        # N of first residue
                        w2_start = np.array(w2_coords[1])
                        # C of third residue
                        w2_end = np.array(w2_coords[11])
                        w2_vec = w2_end - w2_start
                        w2_mag = np.linalg.norm(w2_vec)

                        # $cos(\theta) = \frac{ a \cdot b }{ \| a \| \| b \| }$
                        vector_cosine = np.dot(w1_vec, w2_vec) / w1_mag / w2_mag

                        if vector_cosine < min_cosine:
                            skip = True
                        else:
                            row['cosine'] = vector_cosine

                    # NEIGHBORS: insert ALLOWED
                    if neighbors and not skip:
                        yield row
                # If mp_dist was not met or clash, row is disallowed
                # NON-NEIGHBORS: insert DISALLOWED
                if not neighbors and (not it[0] or skip):
                    yield row
                it.iternext()
            logging.debug("Finished %s %s", w1, w2)

        for windowid in range(1, len(window_indices)):
            window_index = window_indices[windowid]
            window_sql = self.create_modeldist_tables(pdbid=pdbid,
                                                      windowid=windowid,
                                                      windowindex_list=window_indices)
            if window_sql is None: continue
            with shared.write_conn(window_sql['db_name']) as conn:
                cursor = conn.cursor()
                for prev_window_index, sql in window_sql['window_dict'].iteritems():
                    row_gen = calculate_modeldist(prev_window_index, window_index)
                    try:
                        first = next(row_gen)
                    except StopIteration:
                        raise LoadModelScoresError("No allowed pairs for %s %s", window_index, prev_window_index)
                    else:
                        cursor.execute(sql['insert'], first)
                        cursor.executemany(sql['insert'], row_gen)
                    cursor.execute(sql['index'])

    def create_modeldist_tables(self, pdbid, windowid, windowindex_list):
        """
        Create database and tables
        """
        db_name = self.modeldist_sql_data['db_name_fmt'].format(pdbid=pdbid, windowid=windowid)
        if not shared.missing(db_name):
            return
        logging.debug("Creating new DB for %s", windowid)
        allowed_id = windowid - 1
        allowed_sql = self.modeldist_sql(cur_window=windowid,
                                        prev_window=allowed_id,
                                        mode="allowed",
                                        **self.modeldist_sql_data)
        allowed_schema = allowed_sql.pop('schema')
        window_dict = {windowindex_list[allowed_id]: allowed_sql}
        with shared.new_conn(db_name) as conn:
            cursor = conn.cursor()
            cursor.execute(allowed_schema)
            for prev_window in range(windowid - 1):
                disallowed_sql = self.modeldist_sql(cur_window=windowid,
                                                   prev_window=prev_window,
                                                   mode="disallowed",
                                                   **self.modeldist_sql_data)
                cursor.execute(disallowed_sql.pop('schema'))
                window_dict[windowindex_list[prev_window]] = disallowed_sql
        return dict(db_name=db_name,
                    window_dict=window_dict)

    def get_coords(self, row, parser, chain, length=9):
        backbone_atoms = self.params['backbone_atoms']
        cb = self.params['cb']
        stored_atoms = self.params['stored_atoms']
        modelid = row['modelid']
        structure = parser.get_structure(modelid, shared.strip_h(row['path']))
        chainstruc = structure[0][chain]
        # Don't enforce length check if length is None
        if length is not None:
            if len(chainstruc) != length:
                raise LoadModelScoresError("%s has %s residues, expected %s" % (modelid, len(chainstruc), length))
        coordinates = list()
        for x, residue in enumerate(chainstruc):
            het, resi, icode = residue.id
            resname = residue.resname
            if het.strip(): continue
            if icode.strip():
                raise LoadModelScoresError("Can't handle PDB insertion code")
            residue_coordinates = list()
            # Load atoms
            atom_coords = {name: residue[name].coord
                           for name in backbone_atoms}
            if resname == "GLY":
                # Generate virtual CB for GLY
                vectors = {name: residue[name].get_vector()
                           for name in backbone_atoms}
                # Center
                n = vectors['N'] - vectors['CA']
                c = vectors['C'] - vectors['CA']
                # Find rotation axis around N for CA-C vector
                rot = PDB.rotaxis(-np.pi * 120.0 / 180.0, c)
                # Apply rotation to CA-N vector
                cb_at_origin = n.left_multiply(rot)
                # Put on top of CA atom
                CB = cb_at_origin + vectors['CA']
                atom_coords[cb] = CB.get_array()
            else:
                atom_coords[cb] = residue[cb].coord

            for name in stored_atoms:
                coords = [float(c) for c in atom_coords[name]]
                residue_coordinates.append(coords)

            assert len(residue_coordinates) == len(stored_atoms)
            coordinates.extend(residue_coordinates)

        if not coordinates:
            raise LoadModelScoresError("No coordinates")
        coordinates = json.dumps(coordinates)

        row = dict(modelid=modelid, di=row['di'], coordinates=coordinates)
        return pd.Series(row)

    @staticmethod
    def scale_scores(df, iv_list):
        scaled_df = pd.DataFrame()
        for iv_row in iv_list:
            iv = iv_row['columns']
            iv_series = df.loc[:, iv]
            iv_scaled = (iv_series - iv_series.mean()) / iv_series.std()
            if not iv_row['ascending']:
                iv_scaled = iv_scaled * -1
            scaled_df.loc[:, iv] = iv_scaled
        if "dfire" in scaled_df.columns and "itscore" in scaled_df.columns:
            df.loc[:, 'di'] = scaled_df[['dfire', 'itscore']].sum(axis=1)
        return df

    def prepare_sql(self, complexdb):

        sql_kwargs = dict(
            model_tablename="allmodel",
            model_choose_tablename="model",
            model_choose_columns=["modelid", "di", "coordinates"],
        )

        model_choose_schema = """
        CREATE TABLE IF NOT EXISTS {model_choose_tablename}
        (
        {model_choose_columns[0]} INTEGER PRIMARY KEY NOT NULL,
        {model_choose_columns[1]} REAL NOT NULL,
        {model_choose_columns[2]} TEXT NOT NULL,
        FOREIGN KEY({model_choose_columns[0]}) REFERENCES {model_tablename}({model_choose_columns[0]})
        )""".format(**sql_kwargs)

        model_choose_insert = shared.create_insert_statement(sql_kwargs['model_choose_tablename'], sql_kwargs['model_choose_columns'])

        n_model_q_fmt = """SELECT count(*) AS ct
        FROM {model_choose_tablename}
        JOIN {model_tablename} USING({model_choose_columns[0]})
        JOIN fragment f USING(fragmentindex, windowindex)
        WHERE windowindex={{windowindex}}
        """.format(**sql_kwargs)

        window_data_q_fmt = """SELECT {model_choose_columns[0]}, modelindex, fragmentindex,
        windowindex, m.dfire as dfire, itscore
        FROM {model_tablename} m JOIN fragment f USING(fragmentindex, windowindex)
        WHERE windowindex={{windowindex}} AND m.dfire IS NOT NULL
        """.format(**sql_kwargs)

        n_paths_q = "SELECT count(*) FROM {model_tablename}".format(**sql_kwargs)

        with shared.write_conn(complexdb) as conn:
            curs = conn.cursor()
            curs.execute(model_choose_schema)

        return dict(
                    model_insert=model_choose_insert,
                    window_data_q_fmt=window_data_q_fmt,
                    n_model_q_fmt=n_model_q_fmt,
                    n_paths=n_paths_q,
                   )

    _p = dict(
        backbone_atoms=['N', 'CA', 'C'],
        cb='CB',
        bb_threshold=3.0,
        max_clash=0,
        sticky_max_all=6.0,
        sticky_max_any=10.0,
        min_cosine=0.1,
        min_ifourdist=5.2,
        max_ifourdist=13.6,
        min_isixdist=6.5,
        max_isixdist=18.5,
        min_itwelvedist=3.8,
        res_overlap=3,
    )

    modeldist_sql_data = dict(
        tablename_fmt="modeldist{prev_window}{cur_window}",
        allowed_columns=["modela", "modelb", "ifourdist", "mpdist", "cosine", "edgescore"],
        allowed_fmt="""CREATE TABLE {tablename}
        (
            modela INTEGER NOT NULL,
            modelb INTEGER NOT NULL,
            ifourdist REAL NOT NULL,
            mpdist REAL NOT NULL,
            cosine REAL NOT NULL,
            edgescore REAL NOT NULL,
            PRIMARY KEY (modela, modelb)
        );""",
        disallowed_columns=["modela", "modelb", "mpdist"],
        disallowed_fmt="""CREATE TABLE {tablename}
        (
            modela INTEGER NOT NULL,
            modelb INTEGER NOT NULL,
            mpdist REAL NOT NULL,
            PRIMARY KEY (modela, modelb)
        );""",
        index_fmt="CREATE INDEX IF NOT EXISTS aindex{prev_window} ON {tablename} (modela)",
    )

    @staticmethod
    def modeldist_sql(prev_window, cur_window, mode,
                      tablename_fmt, allowed_columns, allowed_fmt,
                      disallowed_columns, disallowed_fmt,
                      index_fmt, **kwargs):
        mode = mode[0].lower()
        tablename = tablename_fmt.format(prev_window=prev_window, cur_window=cur_window)
        schema_fmt = dict(a=allowed_fmt,
                          d=disallowed_fmt)[mode]
        schema = schema_fmt.format(tablename=tablename)
        columns = dict(a=allowed_columns,
                       d=disallowed_columns)[mode]
        insert = shared.create_insert_statement(tablename=tablename,
                                              columns=columns)
        index = index_fmt.format(tablename=tablename, prev_window=prev_window)
        return dict(schema=schema,
                    insert=insert,
                    index=index)

    @property
    def params(self):
        "Compute parameters"
        p = self._p
        p['stored_atoms'] = p['backbone_atoms'] + [p['cb']]
        p['ca_index'] = p['stored_atoms'].index("CA")
        p['n_atoms'] = len(p['stored_atoms'])
        p['atom_overlap'] = p['res_overlap'] * p['n_atoms']
        return p

    @classmethod
    def parse_args(cls, module_args=None):
        a = argparse.ArgumentParser()
        a.add_argument("-d", "--working_dir", metavar="DIR",
                       help="Location of files.")
        a.add_argument("-b", "--complexdb",
                       help="Path to complex database")
        a.add_argument("-p", "--complexname",
                       help="PDB ID (name of subdirectory)")
        a.add_argument("-r", "--receptor_chain",
                       help="Receptor chain")
        a.add_argument("-l", "--ligand_chain",
                       help="Ligand chain")

        args = a.parse_args(module_args)
        return cls(**vars(args))

if __name__ == "__main__":
    LoadModelScores.parse_args()
