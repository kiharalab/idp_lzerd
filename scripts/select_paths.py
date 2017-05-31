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

import argparse
import inspect
import logging
import os

from Bio import PDB
import numpy as np
import pandas as pd
pd.set_option("max_colwidth", 100)

from combine_receptor import CombineChain
import shared

script_dir = os.path.abspath(os.path.dirname(inspect.getfile(inspect.currentframe())))

logging.basicConfig(level=logging.DEBUG)


b_weight = 0.3
neco_weights = [0.5, 0.1, 0.3, 0.1]
neco_scores = (("nodescore", True),
               ("edgescores", True),
               ("clustersize", False),
               ("occupancyscore", False))


class SelectPathsError(RuntimeError):
    """Exception for class SelectPaths"""
    def __init__(self, message):
        super(SelectPathsError, self).__init__(message)


class SelectPaths(object):

    filelist = "filelist.csv"
    default_count = 100

    def __init__(self, complexname, receptor_chain, ligand_chain, nwindows, directory=None, count=None, **kwargs):
        """CONSTRUCTOR"""

        pdb_kwargs = dict(complexname=complexname,
                          receptor_chain=receptor_chain.upper(),
                          ligand_chain=ligand_chain.upper(),
                          nwindows=nwindows)
        pdbwindowid = "{complexname}{receptor_chain}{ligand_chain}{nwindows}".format(**pdb_kwargs)

        if directory is None:
            directory = script_dir

        self.working_dir = os.path.abspath(directory)

        dirname = pdbwindowid.lower()
        local_dir = os.path.join(self.working_dir, dirname)

        if count is None:
            count = self.default_count

        self.select_paths(dest=local_dir, ct=count, **pdb_kwargs)

    def select_paths(self, complexname, receptor_chain, ligand_chain, nwindows, ct, dest=None):

        pdb_kwargs = dict(complexname=complexname,
                          receptor_chain=receptor_chain,
                          ligand_chain=ligand_chain,
                          nwindows=nwindows)
        pdbid = "{complexname}{receptor_chain}{ligand_chain}".format(**pdb_kwargs)
        pdbwindowid = "{0}{nwindows}".format(pdbid, **pdb_kwargs)

        # Create top directory for pdbid
        # Equivalent to directory argument to constructor
        if dest is None:
            # Default dest is pdbid and window number
            dest = os.path.join(self.working_dir, pdbwindowid)
        else:
            # Place non-absolute dest relative to model db filedir
            if not os.path.isabs(dest):
                dest = os.path.join(self.working_dir, dest)

        # charmm balks at mixed case
        dest = dest.lower()
        shared.mkdir_p(dest)

        path_db = "path_{0}_all.db".format(pdbid)
        path_db = os.path.join(self.working_dir, path_db)
        windows = ["window%s" % x for x in range(nwindows)]
        center_q = """SELECT
        pathsid, nodescore, edgescores, clustersize,
        {windows}
        FROM clusters{nwindows}
        JOIN paths{nwindows} USING (pathsid)
        WHERE is_medoid=1
        """.format(nwindows=nwindows,
                   windows=", ".join(windows))

        center_df = shared.db_to_pandas(center_q, path_db)

        occupancy_csv = "{0}_receptor_occupancy.csv".format(pdbwindowid)
        occupancy_file = os.path.join(self.working_dir, occupancy_csv)
        if shared.missing(occupancy_file):
            logging.warning("%s missing", occupancy_file)
            raise SelectPathsError("No occupancy score")
        # Load occupancy score
        occ_data = pd.read_csv(occupancy_file)
        # Combine occupancy score and other scores
        occ_data.rename(columns=dict(pathid="pathsid"), inplace=True)
        center_df = center_df.merge(occ_data, how="left")
        missing = center_df[center_df.isnull().any(axis=1)]
        if not missing.empty:
            print missing
            raise SelectPathsError("Null scores")

        for x, (scorename, ascending) in enumerate(neco_scores): 
            multiplier = 1
            if not ascending:
                multiplier = -1

            if any(pd.isnull(center_df[scorename])):
                logging.error("%s %s", pdbid, scorename)
                raise SelectPathsError("Null values")
            # compute Z-scores
            center_df[scorename + "z"] = self.zscore(center_df[scorename] * multiplier)
        center_df["best_score"] = center_df.apply(lambda x: min(x[s + "z"] for s, __ in neco_scores), axis=1)

        # compute weighted score
        notb_weight = 1 - b_weight
        notb_scores = [(wght, scrnm) for wght, (scrnm, __) in zip(neco_weights, neco_scores)]
        def score_row(r):
            return b_weight * r['best_score'] + notb_weight * sum(wght * r[scrnm + "z"] for wght, scrnm in notb_scores)
        center_df['weighted_score'] = center_df.apply(score_row, axis=1) 
        #print center_df.head()
        # take top n
        sorted = center_df.sort_values('weighted_score')
        top_n = sorted.head(ct)
        top_n[['pathsid', 'nodescorez', 'edgescoresz', 'clustersizez', 'occupancyscorez', 'best_score', 'weighted_score']].to_csv(os.path.join(dest, "path_scores.csv"), index=False)
        paths = top_n.loc[:, ['pathsid'] + windows]

        model_db_file = "scores_{0}.db".format(pdbid)
        model_db_file = os.path.join(self.working_dir, model_db_file)

        self.combine_paths(paths=paths, model_db_file=model_db_file, dest=dest, **pdb_kwargs)

    def combine_paths(self, paths, complexname, receptor_chain, ligand_chain, nwindows, model_db_file, dest):
        """
        Prepare chosen paths for CHARMM relaxation.
        """

        top_dir = dest

        paths = shared.add_model_path(paths, model_db_file)

        pdb_kwargs = dict(complexname=complexname,
                          receptor_chain=receptor_chain,
                          ligand_chain=ligand_chain,
                          nwindows=nwindows)

        window_str = "window"
        window_vars = [x for x in paths.columns.values.tolist() if x.startswith(window_str)]
        window_vars.sort(key=lambda x: int(x.replace(window_str, "")))

        # Get start residue number for each window
        window_start_q = "SELECT res_start FROM window ORDER BY windowindex"
        window_starts = shared.db_to_pandas(window_start_q, model_db_file)['res_start'].values


        parser = PDB.PDBParser(QUIET=True)
        get_structure = lambda x: parser.get_structure(os.path.splitext(os.path.basename(x))[0], shared.strip_h(x))
        io = PDB.PDBIO()
        struct_name = "combined"
        model_dir = os.path.join(self.working_dir, complexname)
        receptor_filename = "{receptor_chain[0]}-{complexname}.pdb".format(**pdb_kwargs)
        receptor_file = os.path.join(model_dir, receptor_filename)
        # if receptor chain is longer than 1, undo combine_chains
        if len(receptor_chain) > 1:
            receptor_file = CombineChain.undo(input=receptor_file)
        # load receptor structure
        receptor_structure = get_structure(receptor_file)
        # reuse receptor structure
        model_id = 0
        receptor_model = receptor_structure[model_id]

        def merge_path(s):
            """
            Create subdirectory and combined.pdb for each path
            """
            # Create subdirectories for pathsid
            pathsid = s['pathsid']
            subdir = os.path.join(top_dir, str(pathsid))
            shared.mkdir_p(subdir)
            outfile = os.path.join(subdir, "%s.pdb" % struct_name)
            if not shared.missing(outfile):
                return outfile
            files = [s[w] for w in window_vars]
            structures = [get_structure(f) for f in files]
            chains = [struc[model_id][ligand_chain] for struc in structures]

            for s_start, c in zip(window_starts, chains):
                # Collect all residues (not modifying chain)
                r_list = [r for r in c]
                # Remove and re-number all residues
                for r in r_list:
                    c.detach_child(r.id)
                    cur_id = list(r.id)
                    cur_id[1] += s_start - 1
                    r.id = tuple(cur_id)
                # Re-add residues to empty chain
                for r in r_list:
                    c.add(r)

            starts = [c.child_list[0].id[1] for c in chains]
            ends = [c.child_list[-1].id[1] for c in chains]

            sb = PDB.StructureBuilder.StructureBuilder()
            sb.init_structure(struct_name)
            sb.init_model(model_id)
            sb.init_seg('    ')
            # Create empty ligand chain
            sb.init_chain(ligand_chain)
            new_struct = sb.get_structure()
            # Add receptor chains
            for ch in receptor_chain:
                new_struct[model_id].add(receptor_model[ch])

            new_chain = new_struct[model_id][ligand_chain]

            for x in xrange(min(starts), max(ends) + 1):
                # Retrieve all residues with id 'x'
                residues = [c[x] for c in chains if x in c]
                # Running total of segment IDs
                n_res = len(residues)
                if n_res == 1:
                    # Unpack single item
                    res, = residues
                    new_chain.add(res)
                elif n_res == 2:
                    # Combined gets averaged position of two residues
                    res1, res2 = residues
                    new_res = res1.copy()
                    for atom1 in res1:
                        atomname = atom1.name
                        atom2 = res2[atomname]
                        new_atom = new_res[atomname]
                        coord1 = atom1.coord
                        coord2 = atom2.coord
                        avg_coord = (coord1 + coord2) / 2.0
                        new_atom.set_coord(avg_coord)
                    new_chain.add(new_res)
                else:
                    raise SelectPathsError("%s residues at %s", n_res, x)

            io.set_structure(new_struct)
            io.save(outfile)
            return outfile

        outfiles = paths.apply(merge_path, 1)

        outfiles.to_csv(os.path.join(top_dir, self.filelist), header=False, index=False)

    @staticmethod
    def zscore(s):
        return (s - np.mean(s)) / np.std(s)

    @classmethod
    def commandline(cls, module_args=None):
        desc = """HELPDESCRIPTION"""
        a = argparse.ArgumentParser(description=desc)
        a.add_argument("complexname",
                       help="4-character PDB code")
        a.add_argument("-r", "--receptor_chain", required=True,
                       help="Receptor chain")
        a.add_argument("-l", "--ligand_chain", required=True,
                       help="Ligand chain")
        a.add_argument("-n", "--nwindows", type=int, required=True,
                       help="Number of windows")
        a.add_argument("-d", "--directory", metavar="DIR",
                       help="Location of files.")
        a.add_argument("-c", "--count", type=int,
                       help="Number of paths to choose")

        args = a.parse_args(module_args)
        c = cls(**vars(args))
        return c


if __name__ == "__main__":
    SelectPaths.commandline()
