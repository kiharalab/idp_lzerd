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
import shutil

from Bio import PDB
from Bio import SeqIO
import pandas as pd

script_dir = os.path.abspath(os.path.dirname(inspect.getfile(inspect.currentframe())))

import shared

logging.basicConfig(level=logging.DEBUG)


class MakePdbError(RuntimeError):
    """Exception for class MakePdb"""
    def __init__(self, message):
        super(MakePdbError, self).__init__(message)


class MakePdb(object):

    header_fmt = "HEADER %-73s\n"
    pdb_fmt = "%-6s%5d %4s%1s%3s %1s%4d%1s   %8s%8s%8s%6s%6s          %2s%2s\n"
    pdb_default = ['ATOM', # Record type
                   None, # Atom number
                   " CA ", # Atom name
                   "", # Altloc
                   None, # Residue name
                   None, # Chain
                   None, # Residue number
                   "", # icode
                   None, # X
                   None, # Y
                   None, # Z
                   "1.00", # occupancy
                   "0.00", # temperature factor
                   "C", # element
                   "" # charge
                  ]

    @classmethod
    def run(cls, fragment_file, ligand_chain, ligand_sequence, pdbid):
        base_dir = os.path.dirname(fragment_file)

        ligand_sequence = SeqIO.read(ligand_sequence, "fasta")

        windowdf = pd.DataFrame(shared.create_windows(len(ligand_sequence)))
        pos_list = windowdf['position'].drop_duplicates().tolist()
        windowdf.to_csv(os.path.join(base_dir, "{0}_data.csv".format(pdbid)), index=False)

        print ligand_sequence

        ### Convert Rosetta format to PDB format
        pos_file_dict = dict()
        with open(fragment_file, "r") as ih:
            position = None
            src_pdbid = None
            rows = list()
            for line in ih:
                parts = line.split()
                if line[:8] == "position":
                    position = parts[1]
                    position_path = os.path.join(base_dir, position)
                    position = int(position)
                    keep_position = position in pos_list
                    if keep_position:
                        pos_file_dict[position] = list()
                        if not os.path.isdir(position_path):
                            os.mkdir(position_path)
                        logging.debug("Rosetta to CA for %s", position_path)
                        index = 1
                # Parts is an empty list if line is blank
                elif not parts:
                    if keep_position and rows:
                        assert src_pdbid
                        filepath = os.path.join(position_path, "frag_%03d.pdb" % index)
                        pos_file_dict[position].append(filepath)
                        if shared.missing(filepath):
                            with open(filepath, "w") as oh:
                                oh.writelines(rows)
                        rows = list()
                        index += 1
                    src_pdbid = None
                elif keep_position:
                    pdbcode, pdbchain, resi, resn, ss, phi, psi, omega, x, y, z = parts
                    new_pdbid = pdbcode + pdbchain
                    if not rows:
                        src_pdbid = new_pdbid
                        rows.append(cls.header_fmt % src_pdbid)
                        res_id = position
                    fmt_list = list(cls.pdb_default)
                    query_idx = res_id - 1
                    assert query_idx >= 0
                    try:
                        query_resn = ligand_sequence[query_idx]
                    except IndexError:
                        print position, index, res_id
                        raise
                    real_res_id = res_id
                    fmt_list[1] = real_res_id
                    fmt_list[4] = shared.one_to_three[query_resn]
                    fmt_list[5] = ligand_chain
                    fmt_list[6] = real_res_id
                    fmt_list[8] = x
                    fmt_list[9] = y
                    fmt_list[10] = z
                    rows.append(cls.pdb_fmt % tuple(fmt_list))
                    res_id += 1

        all_pos = sorted(pos_file_dict.keys(), key=int)
        last_pos = all_pos[-1]

        # Truncate last pos if necessary
        # 1, 7, 13, 19, 25, ... are starts
        if last_pos % 6 != 1:
            parser = PDB.PDBParser(QUIET=True)
            io = PDB.PDBIO()
            # Get computed position from database
            new_start = windowdf[windowdf['position'] == last_pos]['res_start'].tolist()[0]
            assert new_start % 6 == 1
            last_pos_dir = os.path.dirname(pos_file_dict[last_pos][0])
            new_dir = os.path.join(os.path.dirname(os.path.normpath(last_pos_dir)), "{0:.0f}".format(new_start))
            logging.debug("Changing position %s to start at %s", last_pos, new_start)
            shared.mkdir_p(new_dir)
            # ADD NEW DIR TO DICT
            pos_file_dict[new_start] = list()
            residue_remove_slice = slice(new_start - last_pos)
            for fn in pos_file_dict.pop(last_pos):
                structure = parser.get_structure("fragment", fn)
                if len(structure.child_list) != 1:
                    raise MakePdbError("More than one model in %s" % fn)
                model = structure.child_list[0]
                if len(model.child_list) != 1:
                    raise MakePdbError("More than one chain in %s" % fn)
                chain = model[ligand_chain]
                for del_res in chain.get_list()[residue_remove_slice]:
                    chain.detach_child(del_res.id)
                basename = os.path.basename(fn)
                outfile = os.path.join(new_dir, basename)
                io.set_structure(structure)
                io.save(outfile)
                pos_file_dict[new_start].append(outfile)
            shutil.rmtree(last_pos_dir)

    @classmethod
    def parse_args(cls, module_args=None):
        a = argparse.ArgumentParser()
        a.add_argument("-p", "--pdbid",
                       help="PDB ID")
        a.add_argument("-l", "--ligand_chain",
                       help="Ligand chain")
        a.add_argument("-s", "--ligand_sequence",
                       help="FASTA file containing sequence of ligand")
        a.add_argument("-f", "--fragment_file",
                       help="Output of Rosetta fragment picker")

        args = a.parse_args(module_args)
        return cls.run(**vars(args))

if __name__ == "__main__":
    MakePdb.parse_args()

