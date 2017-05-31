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
import subprocess

script_dir = os.path.abspath(os.path.dirname(inspect.getfile(inspect.currentframe())))

import shared

logging.basicConfig(level=logging.DEBUG)


class RunRosettaError(RuntimeError):
    """Exception for class RunRosetta"""
    def __init__(self, message):
        super(RunRosettaError, self).__init__(message)


class RunRosetta(object):

    #nr = "/apps/blast+/databases/nr"
    blastcmdfmt = "{blastpgp_exe} -t 1 -i {id}.fasta -F F -j2 -o {id}.blast -d {nr_path} -v10000 -b10000 -K1000 -h0.0009 -e0.0009 -C {id}.chk -Q {id}.pssm"
    convert_blast = os.path.join(script_dir, "parse.pl")
    ss_methods = ("psipred", "porter", "jpred", "sspro")

    default_nfrag = 30

    @classmethod
    def run(cls, complexname, ligand_chain, 
            ligand_sequence, 
            psipred_path, porter_path, jpred_path, sspro_path,
            directory=None, nfrag=None, **kwargs):

        config = shared.load_config()
        make_fragments_pl = os.path.join(config['rosetta_path'], "tools/fragment_tools/make_fragments.pl")
        fragment_picker_exe = os.path.join(config['rosetta_path'], "main/source/bin/fragment_picker.linuxgccrelease")

        #complexname = complexname[:4]

        if directory is None:
            directory = os.path.join(script_dir, "quota{0}".format(complexname))
        if nfrag is None:
            nfrag = cls.default_nfrag

        directory = os.path.abspath(directory)
        logging.info("DIRECTORY: %s", directory)

        # Check, prepare, run Rosetta
        pdbid = "{0}{1}".format(complexname, ligand_chain)
        output_dir = os.path.join(directory, "output_files")
        fragment_name = "{0}.{1}.9mers".format(pdbid, nfrag)
        fragment_file = os.path.join(output_dir, fragment_name)
        score_name = "{0}.fsc.{1}.9mers".format(pdbid, nfrag)
        score_file = os.path.join(output_dir, score_name)
        input_dir = os.path.join(directory, "input_files")
        path_id = os.path.join(input_dir, pdbid)
        fastain = path_id + ".fasta"
        if shared.missing(fragment_file) or shared.missing(score_file):
            flag_kwargs = dict(pdbid=pdbid, nfrag=nfrag, rosetta_path=config['rosetta_path'])
            template_dir = os.path.join(script_dir, "rosetta_templates")
            # Create Rosetta tree
            shared.mkdir_p(output_dir)
            shared.mkdir_p(input_dir)
            # Check if ss files exist
            for method in cls.ss_methods:
                method_key = "{0}_path".format(method)
                f = locals()[method_key]
                if shared.missing(f):
                    raise RunRosettaError("File %s not found" % f)
                else:
                    flag_kwargs[method_key] = f
            native_line = ""
            protocol_type = ""
            flag_kwargs['native_line'] = native_line
            flag_kwargs['protocol_type'] = protocol_type
            # Copy fasta file
            try:
                shutil.copy(ligand_sequence, fastain)
            except shutil.Error:
                # same file error
                pass
            with shared.CHDIR(input_dir):
                checkpoint_file = "{0}.checkpoint".format(pdbid)
                if shared.missing(checkpoint_file):
                    chk_file = "{0}.chk".format(pdbid)
                    if shared.missing(chk_file):
                        # Run blast
                        cmd = cls.blastcmdfmt.format(id=path_id, **config)
                        cmd = cmd.split()
                        subprocess.check_call(cmd)
                    # Convert to Rosetta checkpoint format
                    #subprocess.check_call([cls.convert_blast, pdbid])
                    subprocess.check_call([cls.convert_blast, make_fragments_pl, fastain, chk_file])
            # Copy quota sizes
            shutil.copy(os.path.join(template_dir, "quota.def"), input_dir)
            # Copy score weights
            weights_name = "quota-protocol.wghts"
            shutil.copy(os.path.join(template_dir, weights_name), input_dir)
            # Create homolog file
            homolog_file = os.path.join(input_dir, "{0}.homolog_vall".format(pdbid))
            with open(homolog_file, "w") as oh:
                oh.write("{0}\n".format(pdbid))
            # Create flag file
            with open(os.path.join(template_dir, "quota-protocol.flags.template")) as ih:
                flags_template = ih.read()
            with open(os.path.join(directory, "quota-protocol.flags"), "w") as oh:
                oh.write(flags_template.format(**flag_kwargs))

            # Run rosetta
            with shared.CHDIR(directory):
                # XXX ask Steve why I need to do this now
                #cmd = "source /usr/local/bio/Modules/default/init/bash; module load rosetta; fragment_picker.linuxgccrelease @quota-protocol.flags"
                #cmd = "module load rosetta; fragment_picker.linuxgccrelease @quota-protocol.flags"
                #bash_cmd = '/bin/bash -c "{0}"'.format(cmd)
                cmd = [fragment_picker_exe, "@quota-protocol.flags"]
                with open("{0}.out".format(pdbid), "w") as oh:
                    #proc = subprocess.Popen(bash_cmd,
                                            #shell=True,
                                            #stdout=oh,
                                            #stderr=oh)
                    proc = subprocess.Popen(cmd, stdout=oh, stderr=oh)
                    returncode = proc.wait()
                    if returncode:
                        raise RunRosettaError("Rosetta exited %s" % returncode)

            # Check if it's actually done
            if shared.missing(fragment_file) or shared.missing(score_file):
                raise RunRosettaError("Rosetta did not finish but exited 0")

    @classmethod
    def parse_args(cls, module_args=None):
        a = argparse.ArgumentParser()
        a.add_argument("complexname",
                       help="PDB ID")
        a.add_argument("-l", "--ligand_chain",
                       help="Ligand chain")
        a.add_argument("-s", "--ligand_sequence",
                       help="FASTA file containing sequence of ligand")
        a.add_argument("--psipred_path",
                       help="PSIPRED output")
        a.add_argument("--porter_path",
                       help="Porter output")
        a.add_argument("--jpred_path",
                       help="Jpred output")
        a.add_argument("--sspro_path",
                       help="SSPro output")
        a.add_argument("-d", "--directory",
                       help="Working directory")
        a.add_argument("-n", "--nfrag",
                       help="Number of fragments to pick")

        args = a.parse_args(module_args)
        try:
            return cls.run(**vars(args))
        except TypeError:
            print args
            raise


if __name__ == "__main__":
    RunRosetta.parse_args()
