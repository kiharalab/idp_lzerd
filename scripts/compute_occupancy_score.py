#!/usr/bin/env python

from __future__ import division

import argparse
import collections
import inspect
import logging
import os

from Bio import PDB
import pandas as pd
pd.set_option("max_colwidth", 100)
import seaborn as sns
sns.set(context="paper")

import shared

script_dir = os.path.abspath(os.path.dirname(inspect.getfile(inspect.currentframe())))

logging.basicConfig(level=logging.DEBUG)


class PlotPathsError(RuntimeError):
    """Exception for class PlotPaths"""
    def __init__(self, message):
        super(PlotPathsError, self).__init__(message)


class PlotPaths(object):

    def __init__(self, complexname, receptor_chain, ligand_chain, directory=None, nwindows=None, **kwargs):
        """CONSTRUCTOR"""

        pdb_kwargs = dict(complexname=complexname,
                          receptor_chain=receptor_chain,
                          ligand_chain=ligand_chain)

        if nwindows < 2:
            raise PlotPathsError("nwindows must be 2 or greater")

        if directory is None:
            directory = script_dir

        plot_kwargs = self.find_db_files(directory=directory, nwindows=nwindows, **pdb_kwargs)
        plot_kwargs.update(pdb_kwargs)

        self.run(**plot_kwargs)

    @classmethod
    def run(cls, dbf, model_db_file, nwindows, query_dict, **kwargs):
        """
        Load data and compute occupancy score

        :param dbf: file containing path information
        :type dbf: str
        :param nwindows: number of windows
        :type nwindows: int
        :param query_dict: SQL queries
        :type dict

        :returns: none
        """
        plot_query = query_dict['plot']
        plot_rows = shared.db_to_pandas(plot_query, dbf)
        if plot_rows.empty:
            print plot_query
            raise PlotPathsError("No rows")

        center_rows = plot_rows.loc[plot_rows.loc[:, 'is_medoid'] == 1]

        size_q = query_dict['clustersize']
        size_rows = shared.db_to_pandas(size_q, dbf)
        if size_rows.isnull().values.any():
            logging.error(size_q)
            raise PlotPathsError("Null sumcsize, may need to rerun find")
        center_rows = center_rows.merge(size_rows, on="pathsid")

        combine_kwargs = dict(nwindows=nwindows, dbf=dbf, model_db_file=model_db_file, query_dict=query_dict)
        for key in "complexname", "receptor_chain", "ligand_chain":
            combine_kwargs[key] = kwargs[key]

        ## Put # of paths each receptor residue contacts into bfactor
        cls.count_receptor_contacts(center_rows, **combine_kwargs)

    @classmethod
    def count_receptor_contacts(cls, paths, complexname, receptor_chain, ligand_chain, nwindows, dbf, model_db_file, query_dict):
        """
        Count number of paths contacting each receptor
        """

        wd = os.path.dirname(model_db_file)
        pdb_kwargs = dict(complexname=complexname,
                          receptor_chain=receptor_chain,
                          ligand_chain=ligand_chain,
                          nwindows=nwindows)
        pdbwindowid = "{complexname}{receptor_chain}{ligand_chain}{nwindows}".format(**pdb_kwargs)
        outfile = os.path.join(wd, "{0}_path_contacts.pdb".format(pdbwindowid))
        path_score_file = os.path.join(wd, "{0}_receptor_occupancy.csv".format(pdbwindowid))
        if not shared.missing(outfile) and not shared.missing(path_score_file):
            logging.debug("%s exists", outfile)
            return

        cutoff = 5.0

        residue_fmt = "{chain}_{resname}{resid[1]}"
        def make_key(residue):
            __, __, chainid, residueid = residue.get_full_id()
            return residue_fmt.format(chain=chainid, resname=residue.get_resname(), resid=residueid)

        # Drop window1 (modelid) column
        orig_window_vars = [x for x in paths.columns.values.tolist() if x.startswith("window")]
        for window_var in orig_window_vars:
            paths = paths.drop(window_var, axis=1)
        # Get model filepaths for paths
        filepaths = cls.get_paths(paths[['pathsid']], dbf=dbf, model_db_file=model_db_file, query_dict=query_dict)
        window_vars = [x for x in filepaths.columns.values.tolist() if x.startswith("window")]
        get_files = lambda row: [row[w] for w in window_vars]

        parser = PDB.PDBParser(QUIET=True)
        # Remove hydrogens
        get_structure = lambda x: parser.get_structure(os.path.splitext(os.path.basename(x))[0], shared.strip_h(x))
        modelid = 0

        receptor_contacts = collections.defaultdict(set)
        for x, row in filepaths.iterrows():
            pathsid = row['pathsid']
            path_files = get_files(row)
            for fn in path_files:
                structure = get_structure(fn)
                atoms = [atom
                         for chain in structure[modelid]
                         for residue in chain
                         for atom in residue]
                if not atoms:
                    raise PlotPathsError("No atoms in %s" % fn)
                ns = PDB.NeighborSearch(atoms)
                search = ns.search_all(radius=cutoff, level="R")
                for res1, res2 in search:
                    __, __, c1, r1 = res1.get_full_id()
                    __, __, c2, r2 = res2.get_full_id()
                    # Skip if chains are both ligand or both receptor
                    if (c1 == ligand_chain) == (c2 == ligand_chain):
                        continue
                    if c1 in receptor_chain:
                        key = make_key(res1)
                    elif c2 in receptor_chain:
                        key = make_key(res2)
                    else:
                        raise PlotPathsError("Neither %s nor %s is receptor" % (c1, c2))
                    receptor_contacts[key].add(pathsid)
        # Convert from defaultdict to normal dict
        receptor_contacts = dict(receptor_contacts)

        # Count paths contacting each receptor residue
        emptyset = set()
        # Chains have been combined
        r_ch = receptor_chain[0]
        # Deliberately using last structure from loop
        for residue in structure[modelid][r_ch]:
            key = make_key(residue)
            mypaths = receptor_contacts.get(key, emptyset)
            count = len(mypaths)
            for atom in residue:
                atom.set_bfactor(count)

        # Write out structure with b-factor
        #structure[modelid].detach_child(ligand_chain)
        #io = PDB.PDBIO()
        #io.set_structure(structure)
        #io.save(outfile)

        # Count receptor contacts for each path
        path_score_dict = collections.defaultdict(int)
        for contacts in receptor_contacts.itervalues():
            n_contacts = len(contacts)
            for pathid in contacts:
                path_score_dict[pathid] += n_contacts
        path_score_df = pd.DataFrame(path_score_dict.items(), columns=["pathid", "occupancyscore"])
        path_score_df.to_csv(path_score_file, index=False)

    @classmethod
    def get_paths(cls, paths, dbf, model_db_file, query_dict):
        paths_query = query_dict['paths']
        allpaths = shared.db_to_pandas(paths_query, dbf)
        # Inner join only keeps path data for given paths
        paths = paths.merge(allpaths, on="pathsid")
        paths = shared.add_model_path(paths, model_db_file=model_db_file)
        return paths

    @staticmethod
    def find_db_files(complexname, receptor_chain, ligand_chain, directory, nwindows=None):

        pdbid = "{complexname}{receptor_chain}{ligand_chain}".format(complexname=complexname,
                          receptor_chain=receptor_chain,
                          ligand_chain=ligand_chain)
        db_filename = "path_{0}_all.db".format(pdbid)
        dbf = os.path.join(directory, db_filename)
        if not os.path.isfile(dbf):
            raise PlotPathsError("File %s not found" % dbf)

        model_db_file = os.path.join(directory, "scores_{0}.db".format(pdbid))

        table_query = "SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'paths%'"
        table_data = shared.db_to_pandas(table_query, dbf)
        table_data['nwindows'] = table_data['name'].apply(lambda x: int(x.replace("paths", "")))
        table_data.sort_values("nwindows", inplace=True)
        # Default do largest number of windows
        if nwindows is None:
            last_row = table_data.tail(1)
            paths_tablename = last_row.name.values[0]
            nwindows = last_row.nwindows.values[0]
        else:
            row = table_data[table_data['nwindows'] == nwindows]
            paths_tablename = row.name.values[0]

        n_kwargs = dict(n=nwindows)
        sql_kwargs = dict(
            paths_tablename=paths_tablename,
            cluster_tablename="clusters{n}".format(**n_kwargs),
            windows=", ".join(["window%s" % x for x in range(nwindows)]),
        )

        select_fmt = """
        SELECT p{nwindows}.pathsid AS pathsid, {sumclusters} AS sumcsize
        FROM paths{nwindows} p{nwindows}
        {pathjoins}
        {clusterjoins}
        WHERE c{nwindows}.is_medoid=1
        """
        sumclusters = " + ".join(["c{n}.clustersize".format(n=n) for n in range(3, nwindows + 1)])
        pathjoins = "\n".join([
            "JOIN paths{n} p{n} USING({windows})".format(
                n=n,
                windows=", ".join(["window{x}".format(x=x)
                                   for x in range(n)]))
                for n in range(3, nwindows)
        ])
        clusterjoins = "\n".join(["JOIN clusters{n} c{n} ON c{n}.pathsid=p{n}.pathsid".format(n=n) for n in range(3, nwindows + 1)])
        if nwindows == 2:
            sumclusters = "c{0}.clustersize".format(nwindows)
            clusterjoins = "JOIN clusters{0} c{0} ON c{0}.pathsid=p{0}.pathsid".format(nwindows)

        cluster_size_q = select_fmt.format(nwindows=nwindows,
                                           sumclusters=sumclusters,
                                           pathjoins=pathjoins,
                                           clusterjoins=clusterjoins)

        plot_q = """SELECT
        pathsid,
        nodescore,
        p.edgescores as edgescores,
        cid,
        clustersize,
        is_medoid
        FROM {cluster_tablename} c
        JOIN {paths_tablename} p USING(pathsid)
        """.format(**sql_kwargs)

        paths_q = "SELECT pathsid, {windows} FROM {paths_tablename}".format(**sql_kwargs)

        return dict(pdbid=pdbid,
                    model_db_file=model_db_file,
                    dbf=dbf,
                    nwindows=nwindows,
                    query_dict=dict(plot=plot_q,
                                    paths=paths_q,
                                    clustersize=cluster_size_q,
                                   ),
                   )

    @classmethod
    def commandline(cls, module_args=None):
        desc = """HELPDESCRIPTION"""
        a = argparse.ArgumentParser(description=desc)
        a.add_argument("complexname",
                       help="PDB code or identifier")
        a.add_argument("-r", "--receptor_chain", required=True,
                       help="Receptor chain(s)")
        a.add_argument("-l", "--ligand_chain", required=True,
                       help="Ligand chain")
        a.add_argument("-n", "--nwindows", type=int, required=True,
                       help="Number of windows.")
        a.add_argument("-d", "--directory", metavar="DIR",
                       help="Location of files.")

        args = a.parse_args(module_args)
        return cls(**vars(args))


if __name__ == "__main__":
    PlotPaths.commandline()
