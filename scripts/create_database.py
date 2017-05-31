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
import glob
import inspect
import os

import shared
from shared import pd

script_dir = os.path.abspath(os.path.dirname(inspect.getfile(inspect.currentframe())))

window_schema = """
CREATE TABLE IF NOT EXISTS window
    (
        windowindex INTEGER PRIMARY KEY NOT NULL,
        window_wd TEXT NOT NULL,
        res_start INTEGER NOT NULL,
        res_end INTEGER NOT NULL,
        position INTEGER NOT NULL,
        timestamp DATE DEFAULT (datetime('now', 'localtime'))
    );
"""

fragment_schema = """
CREATE TABLE IF NOT EXISTS fragment
    (
        windowindex INTEGER NOT NULL,
        fragmentindex INTEGER NOT NULL,
        timestamp DATE DEFAULT (datetime('now', 'localtime')),
        PRIMARY KEY(windowindex, fragmentindex)
    );
"""


class CreateDatabaseError(RuntimeError):
    """Exception for class LoadModelScores"""
    def __init__(self, message):
        super(CreateDatabaseError, self).__init__(message)


def main(directory, pdbid, r_ch, l_ch, input):

    wd = os.path.abspath(directory)
    subdir = os.path.join(wd, pdbid)
    complexid = "{0}{1}{2}".format(pdbid, r_ch, l_ch)
    complexdb = os.path.join(wd, "scores_{0}.db".format(complexid))

    if not shared.missing(complexdb):
        return

    # Initialize window data
    window_data = pd.read_csv(input)
    windows = list()
    fragments = list()
    # Windows are in 2nd-level subdirectories
    for window_dir in glob.iglob(os.path.join(subdir, "*")):
        # Skip 2nd-level files
        if not os.path.isdir(window_dir):
            continue
        subdir, windowindex = os.path.split(window_dir)
        windowindex = windowindex.lower().replace(pdbid.lower(), "")
        try:
            windowindex = int(windowindex)
        except Exception:
            raise CreateDatabaseError("Expected window directory format $PDBID$WINDOWINDEX (e.g. 1bfg1)")
        window_row = dict(windowindex=windowindex, window_wd=window_dir)
        windows.append(window_row)
        # Fragments are in 3rd-level subdirectories
        for fragment_dir in glob.iglob(os.path.join(window_dir, "*")):
            # Skip 3rd-level files
            if not os.path.isdir(fragment_dir):
                continue
            window_dir, fragmentindex = os.path.split(fragment_dir)
            fragment_row = dict(windowindex=windowindex, fragmentindex=fragmentindex)
            fragments.append(fragment_row)

    window_df = pd.merge(window_data, pd.DataFrame(windows), on="windowindex")

    # Create fragment database
    with shared.new_conn(complexdb) as windowconn:
        cursor = windowconn.cursor()
        # Insert windows into database
        cursor.execute(window_schema)
        w_insert = shared.create_insert_statement("window", window_df.columns)
        cursor.executemany(w_insert, window_df.to_dict("records"))
        # Insert fragments into database
        cursor.execute(fragment_schema)
        insert = shared.create_insert_statement("fragment", ["windowindex", "fragmentindex"])
        cursor.executemany(insert, fragments)

def parse_args(module_args=None):
    a = argparse.ArgumentParser()
    a.add_argument("-d", "--directory",
                   help="Working directory")
    a.add_argument("-p", "--pdbid",
                   help="PDB ID (name of subdirectory)")
    a.add_argument("-r", "--r_ch",
                   help="Receptor chain")
    a.add_argument("-l", "--l_ch",
                   help="Ligand chain")
    a.add_argument("-i", "--input",
                   help="Window data in CSV format")

    args = a.parse_args(module_args)
    return main(**vars(args))

if __name__ == "__main__":
    parse_args()
