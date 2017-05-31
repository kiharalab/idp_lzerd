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
import inspect
import itertools
import json
import logging
import os
import shutil
import sqlite3
import sys
import subprocess

script_dir = os.path.abspath(os.path.dirname(inspect.getfile(inspect.currentframe())))
# Root directory of package
ROOTDIR = os.path.normpath(os.path.join(script_dir, os.pardir))

import numpy as np

import shared


class ClusterLRMSD(object):

    def __init__(self, center_coords, n_res, ligand=None):
        """
        :param center_ids: Files making up centers
        :type center_ids: list of 2-tuples
        """

        if ligand is None:
            ligand = 'B'

        self.ligand = ligand

        self.m = n_res

        # Precompute multiplier for RMSD calculation
        self.factor = np.power(self.m, -0.5)
        # Allocate matrix
        coord_matrix = np.empty((len(center_coords), self.m * 3), np.float32)
        # Fill matrix
        for x, (__, coords) in enumerate(center_coords):
            coord_matrix[x,] = list(itertools.chain(*[coords[y + 1] for y in xrange(0, len(coords) - 1, 4)]))

        self.center_ids = [pathsid for pathsid, __ in center_coords]
        self.coord_matrix = coord_matrix

    def add_center(self, pathsid, coords):
        self.coord_matrix = np.vstack((self.coord_matrix, coords))
        self.center_ids.append(pathsid)

    def find_min(self, coords):
        """
        Find center file with closest distance to target file

        :param targetfiles: Files making up path
        :type targetfiles: list

        :returns: 2-tuple (pathsid, lrmsd)
        """
        # Get target structure
        fixed_coord = np.array(list(itertools.chain(*[coords[y + 1] for y in xrange(0, len(coords) - 1, 4)])))

        # Calculate lrmsd values
        lrmsd_vector = np.linalg.norm(self.coord_matrix - fixed_coord, axis=1) * self.factor
        min_lrmsd = np.amin(lrmsd_vector)
        min_pathsid = self.center_ids[np.argmin(lrmsd_vector)]
        return dict(pathsid=min_pathsid,
                    rmsd=min_lrmsd,
                    coords=fixed_coord)


class ClusterPdbError(RuntimeError):
    "Exception for class ClusterPdb"


class ClusterPdb(object):

    default_wd = os.getcwd()
    default_cutoff = 4.0

    def make_sql(self, complexid, nwindows, limit=None):

        path_select_fmt = '''SELECT pathsid,
        {windows} FROM paths{nwindows}
        '''
        path_count_fmt = "SELECT count(*) FROM paths{nwindows}"
        if limit:
            limit_clause = "\nLIMIT {0}".format(limit)
            path_select_fmt += limit_clause
            path_count_fmt += limit_clause

        cluster_tablename="clusters{nwindows}".format(nwindows=nwindows)
        cluster_columns=["pathsid", "cid", "is_medoid", "clustersize"]
        cluster_kwargs = dict(
            nwindows=nwindows,
            cluster_tablename=cluster_tablename,
            cluster_columns=cluster_columns,
            paths_tablename="paths{}".format(nwindows)
        )

        path_count_sql = path_count_fmt.format(**cluster_kwargs)
        cluster_count_sql = """SELECT count(*) FROM {cluster_tablename}
        """.format(**cluster_kwargs)

        cluster_schemas = """
        CREATE TABLE {cluster_tablename}
        (
        {cluster_columns[0]} INTEGER PRIMARY KEY NOT NULL,
        {cluster_columns[1]} INTEGER NOT NULL,
        {cluster_columns[2]} INTEGER NOT NULL,
        {cluster_columns[3]} INTEGER,
        FOREIGN KEY({cluster_columns[0]}) REFERENCES {paths_tablename}({cluster_columns[0]})
        )""".format(**cluster_kwargs).split(";")

        cluster_insert = shared.create_insert_statement(cluster_tablename, cluster_columns[:-1])

        path_select_sql = path_select_fmt.format(
            windows=", ".join("window{0}".format(x) for x in range(nwindows)),
            **cluster_kwargs)

        logging.debug("\n%s", path_select_sql)

        return dict(
                    path_count=path_count_sql,
                    cluster_count=cluster_count_sql,
                    path_select=path_select_sql,
                    cluster_schemas=cluster_schemas,
                    cluster_insert=cluster_insert)

    def __init__(self, complexid, nwindows, directory=None, limit=None):

        config = shared.load_config()
        self.clust_bin = os.path.join(config['lzerd_path'], "LB3Dclust")

        if directory is None:
            directory = script_dir

        path_db_file = os.path.join(directory, "path_{0}_all.db".format(complexid))
        if shared.missing(path_db_file):
            raise ClusterPdbError("DB file %s not found" % path_db_file)
        model_db_file = os.path.join(directory, "scores_{0}.db".format(complexid))
        if shared.missing(model_db_file):
            raise ClusterPdbError("DB file %s not found" % model_db_file)
        logging.debug("\n%s", model_db_file)

        sql_dict = self.make_sql(complexid=complexid, nwindows=nwindows, limit=limit)

        pconn = sqlite3.connect(path_db_file, isolation_level="EXCLUSIVE")
        pcurs = pconn.cursor()
        # Check done
        try:
            cluster_result = pcurs.execute(sql_dict['cluster_count'])
        except sqlite3.OperationalError:
            cluster_count = 0
        else:
            cluster_count = cluster_result.next()[0]
        n = pcurs.execute(sql_dict['path_count']).next()[0]
        done = (n and (cluster_count == n))
        if not done:
            if cluster_count:
                logging.debug("n paths: %s", n)
                logging.debug("n clusters: %s", cluster_count)
                sys.exit(1)
            path_q = sql_dict['path_select']
            row_gen = pcurs.execute(path_q)
            # Convert result tuples to dict of list
            modelid_dict = {int(row[0]): row[1:] for row in row_gen}
            # Start heuristic clustering
            cluster_gen = self.partial_cluster(modelid_dict=modelid_dict, complexid=complexid, nwindows=nwindows, model_db_file=model_db_file)
            # Create cluster table
            for stmt in sql_dict['cluster_schemas']:
                pcurs.execute(stmt)
            # Insert cluster rows
            insert = sql_dict['cluster_insert']
            # Write to disk in batches
            for path_chunk in itertools.izip_longest(*[iter(cluster_gen)] * 10000, fillvalue=None):
                pcurs.executemany(insert,
                                  (row for row in path_chunk if row is not None))
            pconn.commit()
        else:
            logging.debug("Clustering done.")
        pconn.close()

    def parse_cluster_out(self, out):
        """
        Parse output of clustering program.

        :param out: output of clustering program
        :type out: str

        :returns: tuple of (list of dict, list of str)
        """
        centr_start = "CID="
        child_start = "chil"
        # Indices of columns to keep
        indices_dict = {centr_start: (1, 7),
                        child_start: (2, 6)}
        # Column names
        headers_dict = {centr_start: ("cid", "model"),
                        child_start: ("cid", "model")}
        type_dict = {centr_start: 1,
                     child_start: 0}
        data_list = list()
        err_list = list()
        for row in out.splitlines():
            start = row[:4]
            if start in (centr_start, child_start):
                data_row = dict(zip(headers_dict[start],
                                    (row.split()[x] for x in indices_dict[start])))
                data_row['is_medoid'] = type_dict[start]
                data_list.append(data_row)
            else:
                err_list.append(row + "\n")

        return data_list, err_list

    def partial_cluster(self, modelid_dict, complexid, nwindows, model_db_file, cutoff=None, cleanup=True):
        """
        Cluster a subset of the data and assign the rest of the data to the subsets.

        :param clusters: list
        """
        if cutoff is None:
            cutoff = self.default_cutoff
        default_fraction = 0.1
        max_sample = 10000
        n = len(modelid_dict)
        half = int(n * 0.5)
        fraction = int(n * default_fraction)
        # Use half if n is very small
        if half <= max_sample:
            k = half
        # Use max is n is very large
        elif fraction > max_sample:
            k = max_sample
        # Otherwise use default fraction
        else:
            k = fraction

        # Select k random items
        chosen_keys = np.random.choice(n, k, replace=False)
        # Shift from 0-index to 1-index
        chosen_keys += 1
        chosen_keys = chosen_keys.tolist()
        sample_dict = {pathsid: modelid_dict[pathsid] for pathsid in chosen_keys}
        all_keys = set(range(1, n + 1))
        other_keys = all_keys.difference(chosen_keys)
        other_dict = {pathsid: modelid_dict[pathsid] for pathsid in other_keys}

        model_query = '''
        SELECT modelid, f.windowindex, f.fragmentindex, m.modelindex, coordinates
        FROM model JOIN allmodel m USING(modelid) JOIN fragment f USING(windowindex, fragmentindex)
        '''
        mconn = sqlite3.connect(model_db_file, isolation_level="EXCLUSIVE")
        mcurs = mconn.cursor()
        model_gen = mcurs.execute(model_query)
        headers = ["windowindex", "fragmentindex", "modelindex", "coordinates"]
        model_data = {row[0]: dict(zip(headers, row[1:])) for row in model_gen}

        chain = complexid[-1]
        modelcoord_dict = {pathsid: list(itertools.chain(*[json.loads(model_data[modelid]['coordinates'])
                                                         for modelid in sample_dict[pathsid]]))
                         for pathsid, window_list in sample_dict.iteritems()}

        clusters = self.cluster_models(modelcoord_dict, chain=chain,
                                      groupid=complexid,
                                      cleanup=True)
        max_cluster_id = 0
        center_coords = list()
        center_dict = dict()
        for row in clusters:
            # Yield clusters for database
            yield row
            # Prepare centers
            if row['is_medoid'] == 1:
                new_cluster_id = int(row['cid'])
                if new_cluster_id > max_cluster_id:
                    max_cluster_id = new_cluster_id
                pathsid = row['pathsid']
                # Store rows for later modification
                center_dict[pathsid] = row
                # Load coords
                center_coords.append(
                    [pathsid,
                     list(itertools.chain(*[json.loads(model_data[modelid]['coordinates'])
                                            for modelid in sample_dict[pathsid]]))])

        nres_query = "SELECT count(*), max(res_end) - min(res_start) + 1 FROM window"
        total_windows, nres = mcurs.execute(nres_query).next()
        if nwindows == total_windows:
            # Add 3 for each overlap
            n_res = nres + (total_windows - 1) * 3
        else:
            n_res = nwindows * 9
        lrmsd = ClusterLRMSD(center_coords=center_coords, n_res=n_res, ligand=chain)
        # Find closest cluster center for each non-chosen file
        for pathsid, window_ids in other_dict.iteritems():
            # Load target coords
            target_coords = list(itertools.chain(*[json.loads(model_data[modelid]['coordinates'])
                                                   for modelid in window_ids]))
            min = lrmsd.find_min(coords=target_coords)
            min_rmsd = min['rmsd']
            # Add model to cluster of closest center
            if min_rmsd <= cutoff:
                # Copy center row
                row_dict = center_dict[min['pathsid']].copy()
                # Update new row with current file data
                row_dict['pathsid'] = pathsid
                row_dict['is_medoid'] = 0
                # Add row to clusters
                yield row_dict
            # Distant model becomes new cluster center
            else:
                # Increment cluster id
                new_cluster_id = max_cluster_id + 1
                max_cluster_id = new_cluster_id
                # Create new row
                row_dict = dict(pathsid=pathsid,
                                cid=new_cluster_id,
                                is_medoid=1)
                # Add row to clusters
                yield row_dict
                # Add row to centers
                center_dict[pathsid] = row_dict
                # Add center to coord matrix
                lrmsd.add_center(pathsid=pathsid, coords=min['coords'])

    def cluster_models(self, modelcoord_dict, chain, groupid, cutoff=None, wd=None, cleanup=True):
        """
        Cluster models.

        :param modelrow_dict: paths to cluster keyed by pathsid
        :type modelrow_dict: dict
        :param chain: chain to keep in model
        :param groupid: group identifier
        :type groupid: int
        :param wd: working directory
        :param cleanup: whether to remove merged files
        :type cleanup: bool

        :returns: list of dict
        """
        if cutoff is None:
            cutoff = self.default_cutoff
        if wd is None:
            wd = self.default_wd
        os.chdir(wd)
        ligand_list = os.path.join(wd, "ligand_list_%s.txt" % groupid)
        cluster_err = os.path.join(wd, "clusters_%s_out.txt" % groupid)

        """
        1         2         3         4         5         6         7         8
        12345678901234567890123456789012345678901234567890123456789012345678901234567890
        ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 11.92      A1   N
        ATOM    146  CA  VAL A  25      31.132  16.439  58.160  1.00 11.85      A1   C
        """
        pdb_line = "ATOM  {index:5d}  CA  UNK A{index:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  \n"

        outdir = os.path.join(wd, "{0}merged".format(groupid))
        shared.mkdir_p(outdir)
        file_list = list()
        lig_file_dict = dict()
        for pathsid, coords in modelcoord_dict.iteritems():
            ligandfile = os.path.join(outdir, "path_{0}.pdb".format(pathsid))
            outlines = list()
            for i, (x, y, z) in enumerate(coords):
                outlines.append(pdb_line.format(index=i + 1, x=x, y=y, z=z))
            with open(ligandfile, "w") as oh:
                oh.writelines(outlines)
            lig_file_dict[ligandfile] = pathsid
            file_list.append(ligandfile)
        with open(ligand_list, "w") as oh:
            for fn in file_list:
                oh.write(fn + "\n")

        logging.debug("Clustering...")
        clst_cmd = [self.clust_bin,
                    "-L", ligand_list,
                    "-c", str(cutoff),
                    "-r", "0.1"]
        proc = subprocess.Popen(clst_cmd,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        out, err = proc.communicate()
        ret = proc.returncode
        if ret:
            raise ClusterPdbError("Clustering exited %s" % ret)
        data_rows, err_lines = self.parse_cluster_out(out)
        with open(cluster_err, "w") as eh:
            eh.writelines(err_lines)
        for row in data_rows:
            row['pathsid'] = lig_file_dict[row['model']]

        if cleanup:
            shutil.rmtree(outdir, ignore_errors=True)
            shared.silent_remove(ligand_list)

        return data_rows

    @classmethod
    def commandline(cls, module_args=None):
        desc = """Heuristic clustering of paths"""
        a = argparse.ArgumentParser(description=desc)
        a.add_argument("complexid",
                       help="PDB code + receptor + ligand.")
        a.add_argument("-n", "--nwindows", type=int, required=True,
                       help="Number of windows.")
        a.add_argument("-d", "--directory", metavar="DIR",
                       help="Location of db files.")
        a.add_argument("-v", "--verbose", action="store_true",
                       help="Verbose mode.")
        a.add_argument("-l", "--limit",
                       help="Size limit")

        kwargs = vars(a.parse_args(module_args))
        verbose = kwargs.pop('verbose')
        if verbose:
            print "Setting verbose"
            level = logging.DEBUG
        else:
            level = logging.INFO
        logging.basicConfig(level=level)

        return cls(**kwargs)

if __name__ == "__main__":
    ClusterPdb.commandline()
