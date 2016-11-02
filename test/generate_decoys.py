#!/usr/bin/env python

import errno
import glob
import inspect
import logging
import os
import subprocess
import tempfile

logging.basicConfig(level=logging.DEBUG)

script_dir = os.path.abspath(os.path.dirname(inspect.getfile(inspect.currentframe())))
# Root directory of package
ROOTDIR = os.path.normpath(os.path.join(script_dir, os.pardir))

SUBDIR = "4ah2"
# receptor filename: 4ah2/a-4ah2_h.pdb
RECEPTOR_CHAIN = "a"
# ligand filename(s): 4ah2/4ah2*/*/c-4ah2_h.pdb
LIGAND_CHAIN = "c"


class TestError(RuntimeError):
    """Exception for class Test"""
    def __init__(self, message):
        super(TestError, self).__init__(message)


class TemporaryFile(object):

    def __init__(self, final_path):
        tmpfile_dir = os.path.dirname(final_path)

        self.tmpfile = tempfile.NamedTemporaryFile(dir=tmpfile_dir)
        self.final_path = final_path

    def __getattr__(self, attr):
        return getattr(self.tmpfile, attr)

    def __enter__(self):
        self.tmpfile.__enter__()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is None:
            self.tmpfile.delete = False
            result = self.tmpfile.__exit__(exc_type, exc_val, exc_tb)
            os.rename(self.tmpfile.name, self.final_path)
            # NamedTemporaryFiles are always 0600
            os.chmod(self.final_path, 0o644)
        else:
            result = self.tmpfile.__exit__(exc_type, exc_val, exc_tb)
        return result


def main():

    wd = os.getcwd()
    subdir = os.path.join(wd, SUBDIR)

    # Parse paths
    path_file = os.path.join(ROOTDIR, "PATHS.ini")
    if not os.path.isfile(path_file):
        raise TestError("Could not open PATHS.ini")
    config = dict()
    with open(path_file) as ih:
        for line in ih:
            line = line.strip()
            parts = line.split(":")
            if len(parts) == 2:
                key, val = parts
                key = key.strip()
                val = val.strip()
                config[key] = val
            else:
                if line[0] != "[" or line[-1] != "]":
                    logging.debug("Could not interpret config line %s", line)

    try:
        lzerd_path = config['lzerd_path']
    except KeyError:
        raise TestError("'{0}' did not contain required key 'lzerd_path'".format(path_file))
    else:
        lzerd_filenames = ("charmm_atom_types.in", "charmm_params.in", "uniCHARMM")
        data_files = [os.path.join(lzerd_path, fn) for fn in lzerd_filenames]
        pdbgen_bin = os.path.join(lzerd_path, "PDBGEN")

    if not os.path.isfile(pdbgen_bin):
        logging.debug(pdbgen_bin)
        logging.debug(os.path.isfile(pdbgen_bin))
        raise TestError("Could not find PDBGEN in {0}".format(lzerd_path))

    generate_decoys(subdir=subdir,
                    pdbgen_bin=pdbgen_bin,
                    data_files=data_files)

def generate_decoys(subdir, pdbgen_bin, data_files):
    receptor_file = os.path.join(subdir, "{1}-{0}_h.pdb".format(SUBDIR, RECEPTOR_CHAIN))
    with open(receptor_file) as ih:
        receptor_list = [line for line in ih if line[:4] == "ATOM"]
    pdb_line_fmt = "{:<80}\n"
    ter_line = pdb_line_fmt.format("TER")
    end_line = pdb_line_fmt.format("END")
    receptor_list.append(ter_line)

    for filepath in data_files:
        dest_path = os.path.basename(filepath)
        if not os.path.isfile(dest_path):
            os.symlink(filepath, dest_path)

    errors = list()
    # Fragments are in 3rd-level subdirectories
    for fragment_dir in glob.iglob(os.path.join(subdir, "*", "*")):
        # Skip 3rd-level files
        if not os.path.isdir(fragment_dir):
            continue
        print fragment_dir

        lzerd_output = os.path.join(fragment_dir, "{0}-{1}.cluster4".format(RECEPTOR_CHAIN, LIGAND_CHAIN))

        number = 50000
        outdir = os.path.join(fragment_dir, "decoys")
        try:
            os.mkdir(outdir)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(outdir):
                pass
            else:
                logging.debug("Couldn't create %s", outdir)
                errors.append(fragment_dir)
                continue

        with open(lzerd_output) as ih:
            number = min(sum(1 for line in ih), number)

        n_decoy_models = len(glob.glob(os.path.join(outdir, "model*.pdb")))
        if n_decoy_models >= number:
            logging.info("Model files complete")
            continue

        ligand_file = os.path.join(fragment_dir, "{1}-{0}_h.pdb".format(SUBDIR, LIGAND_CHAIN))

        cmd = [pdbgen_bin, receptor_file, ligand_file,
               lzerd_output, str(number)]
        with open(os.path.join(fragment_dir, "pdbgen.log"), "w") as oh:
            ret = subprocess.call(cmd, stdout=oh, stderr=subprocess.STDOUT) 
            if ret:
                logging.debug("PDBGEN returned %s", ret)
                errors.append(fragment_dir)
                continue

        for x in range(1, number + 1):
            ligand_file = "ligand{}.pdb".format(x)
            model_file = "model{}.pdb".format(x)
            model_file = os.path.join(outdir, model_file)
            # Copy list
            outlist = list(receptor_list)
            with open(ligand_file) as ih:
                outlist.extend(ih)
            outlist.append(ter_line)
            outlist.append(end_line)
            with TemporaryFile(model_file) as oh:
                oh.writelines(outlist)
            os.remove(ligand_file)

    if errors:
        raise TestError("Errors on %s directories" % len(errors))
    else:
        logging.info("Job complete with 0 errors")

if __name__ == "__main__":
    main()
