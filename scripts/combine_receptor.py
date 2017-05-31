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
import math
import os
import string

from Bio import PDB
from Bio.PDB.Chain import Chain

import shared

logging.basicConfig(level=logging.DEBUG)

script_dir = os.path.abspath(os.path.dirname(inspect.getfile(inspect.currentframe())))


class CombineChainError(RuntimeError):
    """Exception for class CombineChain"""
    def __init__(self, message):
        super(CombineChainError, self).__init__(message)


class CombineChain(object):

    allowed_chains = string.ascii_letters + string.digits

    def __init__(self, input, receptor=None, ligand=None):
        """CONSTRUCTOR"""
        if receptor is None:
            receptor = "A"

        if ligand is not None and len(ligand) != 1:
            raise CombineChainError("Ligand can only be one chain: %s" % ligand)

        nchains = len(receptor)
        receptor_chains = [c for c in receptor if c in self.allowed_chains]
        if len(receptor_chains) != nchains:
            logging.info("Changed receptor input %s to %s", receptor, "".join(receptor_chains))

        parser = PDB.PDBParser(QUIET=True)

        structure = parser.get_structure(os.path.splitext(os.path.basename(input))[0], shared.strip_h(input))
        # Remove all other models
        for mdl in structure.get_list()[1:]:
            structure.detach_child(mdl.get_id())

        model = structure[0]
        model_chains = set(c.id for c in model)

        # Remove ligand/receptor from chain set
        error = False
        if ligand is not None:
            try:
                model_chains.remove(ligand)
            except KeyError:
                error = True
        for r in receptor_chains:
            try:
                model_chains.remove(r)
            except KeyError:
                error = True
        if error:
            logging.error("Chains: %s", "".join(c.id for c in model.child_list))
            logging.error("Receptor: %s", "".join(receptor_chains))
            logging.error("Ligand: %s", ligand)
            raise CombineChainError("Specified chains not found in %s" % input)

        # Remove all other chains
        for chainid in model_chains:
            model.detach_child(chainid)

        # Remove heteroatoms
        for chain in model:
            remove = list()
            for res in chain:
                res_id = res.get_id()
                if res_id[0] != " ":
                    remove.append(res_id)
            for res_id in remove:
                chain.detach_child(res_id)
        
        # Renumber receptor chains
        chain_new_old_dict = dict()
        first_chain = None
        prev_chain_end = None
        for chainid in receptor_chains:
            cur_chain = model[chainid]
            if prev_chain_end is None:
                first_chain = cur_chain
            # Add multiples of 100 to make first residue be greater than previous chain end
            else:
                cur_start = cur_chain.child_list[0].get_id()[1]
                # Round offset up to nearest 100
                raw_offset = prev_chain_end - cur_start + 1
                pct_offset = raw_offset * 0.01
                ceil_offset = math.ceil(pct_offset)
                offset = int(ceil_offset * 100)
                # Store offset
                chain_new_old_dict[cur_chain.id] = (cur_start + offset, cur_start)
                # Add offset to all residue numbers
                for res in cur_chain:
                    cur_id = list(res.get_id())
                    cur_id[1] += offset
                    res.id = tuple(cur_id)
                    # Transfer renumbered residue to first chain
                    first_chain.add(res)
                # Remove current chain
                model.detach_child(cur_chain.id)
            # Store resseq of last residue in chain
            prev_chain_end = cur_chain.child_list[-1].get_id()[1]

        pdbcode = os.path.splitext(os.path.basename(input))[0]
        pdbid = "{pdbcode}{rchains}".format(pdbcode=pdbcode,
                                            rchains="".join(receptor_chains))
        if ligand is not None:
            pdbid += ligand
        output_file = "{0}.pdb".format(pdbid)
        io = PDB.PDBIO()
        io.set_structure(structure)
        with open(output_file, "w") as oh:
            for chain_id, old_new in chain_new_old_dict.iteritems():
                line = "REMARK RESIDUE={1[0]} CHAIN={0} START={1[1]}\n".format(chain_id, old_new)
                oh.write(line)
            io.save(oh)
        self.output_file = output_file

    @staticmethod
    def extract_residue_dict(input):
        residue_dict = dict()
        with open(input) as ih:
            for line in ih:
                if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("MODEL"):
                    break
                if line.startswith("REMARK RESIDUE="):
                    parts = line.split()[1:]
                    data = dict(p.split("=") for p in parts)
                    residue_dict[int(data['RESIDUE'])] = data
        return residue_dict

    @classmethod
    def undo(cls, input, outfile=None, ligand_chain=None, residue_dict=None, write=True):
        if residue_dict is None:
            residue_dict = cls.extract_residue_dict(input)
        if not residue_dict:
            raise CombineChainError("No residue dict")
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure(os.path.splitext(os.path.basename(input))[0], input)
        if len(structure.child_list) > 1:
            raise CombineChainError("Input has more than one model: %s" % structure.child_list)
        model = structure.child_list[0]
        chain_id_list = [c.id for c in model.child_list]
        chain_id_set = set(chain_id_list)
        if len(chain_id_list) != len(chain_id_set):
            raise CombineChainError("Non-unique chain IDs")
        if ligand_chain is not None:
            try:
                chain_id_set.remove(ligand_chain)
            except ValueError:
                raise CombineChainError("Ligand chain %s not found in %s", ligand_chain, input)
        if len(chain_id_set) > 1:
            logging.debug(input)
            raise CombineChainError("Input has more than one chain: %s" % model.child_list)
        receptor_chain = chain_id_set.pop()
        chain = model[receptor_chain]
        residue_list = sorted(residue_dict.keys())
        for x, residue in enumerate(residue_list):
            data = residue_dict[residue]
            chainid = data['CHAIN']
            try:
                next_change = residue_list[x + 1]
            except IndexError:
                next_change = None
            new_chain = Chain(chainid)
            model.add(new_chain)
            for res in chain.get_list():
                resseq = res.id[1]
                if (next_change is None or resseq < next_change) and resseq >= residue:
                    chain.detach_child(res.id)
                    new_chain.add(res)
        if write:
            if outfile is None:
                fileparts = list(os.path.splitext(input))
                fileparts.insert(1, "_split")
                outfile = "".join(fileparts)
            if not os.path.isabs(outfile):
                outdir = os.path.dirname(input)
                outfile = os.path.join(outdir, outfile)
            io = PDB.PDBIO()
            io.set_structure(structure)
            io.save(outfile)
            return outfile
        else:
            return structure
        assert False

    @classmethod
    def commandline(cls, module_args=None):
        desc = """HELPDESCRIPTION"""
        a = argparse.ArgumentParser(description=desc)
        a.add_argument("input",
                       help="Input file")
        a.add_argument("-r", "--receptor",
                       help="Receptor chain(s), multiple without spaces e.g. AB")
        a.add_argument("-l", "--ligand",
                       help="Ligand chain")

        args = a.parse_args(module_args)
        c = cls(**vars(args))
        return c


if __name__ == "__main__":
    CombineChain.commandline()
