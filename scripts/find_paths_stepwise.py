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
import logging
import os
import time

import apsw
import numpy as np
import pandas as pd

import shared

script_dir = os.path.abspath(os.path.dirname(inspect.getfile(inspect.currentframe())))

from cluster_heuristic import ClusterPdb

logging.basicConfig(level=logging.DEBUG)

total_paths = 0
skipped_paths = 0


class FindPathsStepwiseError(RuntimeError):
    """Exception for class FindPathsStepwise"""
    def __init__(self, message):
        super(FindPathsStepwiseError, self).__init__(message)


class FindPathsStepwise(object):

    default_batchsize = 10000

    out_db_fmt = os.path.join("{db_dir}", "path_{complexid}_all.db")
    in_tablename_fmt = "modeldist{prev_window}{cur_window}"
    out_tablename_fmt = "paths{nwindows}"
    cluster_tablename_fmt = "clusters{nwindows}"
    model_score_columns = ["di"]
    pair_score_columns = ["edgescore"]
    path_score_columns = ["edgescores", "nodescore"]
    cluster_update_columns = ["clustersize"]
    cluster_id_columns = ["pathsid"]

    def __init__(self, complexid, nwindows, batchsize=None, directory=None):
        """CONSTRUCTOR"""

        if directory is None:
            db_dir = script_dir
        else:
            db_dir = directory
        self.in_db_fmt = os.path.join(db_dir, "%s_modeldist{cur_window}.db" % complexid)
        self.score_db_file = os.path.join(db_dir, "scores_{complexid}.db".format(complexid=complexid))
        self.out_db_file = self.out_db_fmt.format(db_dir=db_dir, complexid=complexid)
        logging.debug(self.out_db_file)

        if batchsize is None:
            batchsize = self.default_batchsize

        make_path_times = list()
        cluster_times = list()
        clustersize_times = list()
        for n in range(2, nwindows + 1):
            ct = self.count_paths(complexid=complexid,
                                  nwindows=n,
                                  db_dir=db_dir)
            n_paths = ct['path_count']
            n_clusters = ct['cluster_count']
            start = time.time()
            if n_paths:
                logging.debug("%s %s has %s rows", complexid, n, n_paths)
            else:
                self.make_paths(complexid=complexid,
                                nwindows=n,
                                batchsize=batchsize)
            pathtime = time.time()
            make_path_times.append(pathtime - start)
            if not n_clusters:
                ClusterPdb(complexid=complexid,
                           nwindows=n,
                           directory=db_dir)
            clustertime = time.time()
            cluster_times.append(clustertime - pathtime)
            # Update cluster sizes
            cluster_tablename = self.cluster_tablename_fmt.format(nwindows=n)
            columns = ["pathsid", "cid", "is_medoid"]
            q = "SELECT {columns} FROM {cluster_tablename}".format(columns=", ".join(columns),
                                                                   cluster_tablename=cluster_tablename)
            # NB shared.db_to_pandas raises factually incorrect error
            with shared.ro_conn(self.out_db_file) as out_conn:
                rows = list(out_conn.cursor().execute(q))
                cluster_rows = pd.DataFrame(rows, columns=columns)
            cluster_sizes = cluster_rows.groupby('cid').size()
            cluster_sizes = cluster_sizes.to_frame("clustersize")
            cluster_rows = cluster_rows.merge(cluster_sizes, left_on="cid", right_index=True)
            center_rows = cluster_rows.loc[cluster_rows.loc[:, 'is_medoid'] == 1]
            cluster_update = shared.create_update_statement(tablename=cluster_tablename,
                                                          columns=self.cluster_update_columns,
                                                          where=self.cluster_id_columns)
            with shared.write_conn(self.out_db_file) as conn:
                conn.cursor().executemany(cluster_update, center_rows.to_dict("records"))
            clustersizetime = time.time()
            clustersize_times.append(clustersizetime - clustertime)
            logging.info("Ending n=%s after %s", n, clustersizetime - start)

        logging.info(make_path_times)
        logging.info(cluster_times)
        logging.info(clustersize_times)
        make_path_total = sum(make_path_times)
        cluster_total = sum(cluster_times)
        clustersize_total = sum(clustersize_times)
        grand_total = sum([make_path_total, cluster_total, clustersize_total])
        logging.info("Make path: %s", make_path_total / grand_total)
        logging.info("Cluster: %s", cluster_total / grand_total)
        logging.info("Clustersize: %s", clustersize_total / grand_total)

    @classmethod
    def count_paths(cls, complexid, nwindows, db_dir):
        out_db_file = cls.out_db_fmt.format(db_dir=db_dir, complexid=complexid)
        out_tablename = cls.out_tablename_fmt.format(nwindows=nwindows)
        cluster_tablename = cls.cluster_tablename_fmt.format(nwindows=nwindows)

        path_count = 0
        cluster_count = 0
        if os.path.isfile(out_db_file):
            # Need write ability to resolve journal
            with shared.write_conn(out_db_file) as out_conn:
                out_conn.cursor().execute("SELECT 1")
            with shared.ro_conn(out_db_file) as out_conn:
                try:
                    path_count = list(out_conn.cursor().execute("SELECT count(*) FROM %s" % out_tablename))[0][0]
                except apsw.SQLError:
                    pass
                try:
                    cluster_count = list(out_conn.cursor().execute("SELECT count(*) FROM %s" % cluster_tablename))[0][0]
                except apsw.SQLError:
                    pass
        return dict(path_count=path_count, cluster_count=cluster_count)

    def make_paths(self, complexid, nwindows, batchsize):
        """Make paths for specific window size"""

        global total_paths
        global skipped_paths

        cur_window = nwindows - 1
        n_edges = cur_window
        prev_window = cur_window - 1
        n_prev_edges = prev_window

        out_db_file = self.out_db_file
        in_db_file = self.in_db_fmt.format(cur_window=cur_window)
        score_db_file = self.score_db_file

        out_tablename = self.out_tablename_fmt.format(nwindows=nwindows)

        out_column_list = ["window%s" % x for x in range(nwindows)]

        out_schema = ["CREATE TABLE {tablename} (".format(tablename=out_tablename)]
        out_schema.append("pathsid INTEGER PRIMARY KEY NOT NULL,")
        out_schema.extend("%s INTEGER NOT NULL," % w_col for w_col in out_column_list)
        out_schema.extend("%s REAL NOT NULL," % col for col in self.path_score_columns)
        out_schema.append("timestamp DATE DEFAULT (datetime('now', 'localtime')));")
        out_schema = "\n".join(out_schema)

        logging.debug("Out schema: %s", out_schema)

        out_column_list.extend(self.path_score_columns)
        insert_sql = shared.create_insert_statement(tablename=out_tablename, columns=out_column_list)
        logging.debug("Out db: %s", out_db_file)
        logging.debug("Insert SQL: %s", insert_sql)

        if nwindows == 2:
            start_db_file = self.in_db_fmt.format(cur_window=1)
            start_table = self.in_tablename_fmt.format(prev_window=0,
                                                       cur_window=1)
            window_columns = ["window%s" % (x) for x in range(2)]
            start_columns = window_columns + self.pair_score_columns
            edgescore_index = start_columns.index("edgescore")
            start_query = """
            SELECT modela AS {windows[0]}, modelb AS {windows[1]}, {score_columns} FROM {tablename}
            """.format(tablename=start_table,
                       windows=window_columns,
                       score_columns=", ".join(self.pair_score_columns))
        elif nwindows == 3:
            start_db_file = self.in_db_fmt.format(cur_window=1)
            start_table = self.in_tablename_fmt.format(prev_window=0,
                                                       cur_window=1)
            window_columns = ["window%s" % (x) for x in range(2)]
            start_columns = window_columns + self.pair_score_columns
            edgescore_index = start_columns.index("edgescore")
            start_query = """
            SELECT modela AS {windows[0]}, modelb AS {windows[1]}, {score_columns} FROM {tablename}
            """.format(tablename=start_table,
                       windows=window_columns,
                       score_columns=", ".join(self.pair_score_columns))
        else:
            start_db_file = out_db_file
            start_table = self.out_tablename_fmt.format(nwindows=cur_window)
            window_columns = ["window%s" % (x) for x in range(cur_window)]
            start_columns = window_columns + self.path_score_columns
            edgescore_index = start_columns.index("edgescores")
            start_query = "SELECT {columns} FROM {tablename}".format(columns=", ".join(start_columns), tablename=start_table)
            cluster = self.cluster_tablename_fmt.format(nwindows=cur_window)
            start_query += " JOIN {clustertable} USING (pathsid) WHERE is_medoid=1".format(clustertable=cluster)

        logging.debug("Starting path db: %s", start_db_file)
        logging.debug("Starting paths query: %s", start_query)
        logging.debug("Edgescore index: %s", edgescore_index)

        logging.debug("Current path db: %s", in_db_file)

        new_window = "window{cur_window}".format(cur_window=cur_window)
        new_edgescore_index = self.pair_score_columns.index("edgescore") + 1
        child_q_fmt = """SELECT modelb AS {new_window}, {pair_score_columns} FROM modeldist{prev_window}{cur_window}
        WHERE modela={{path[{prev_window}]}}""".format(new_window=new_window, cur_window=cur_window, prev_window=prev_window, pair_score_columns=", ".join(self.pair_score_columns))
        for prior_window in range(prev_window):
            child_q_fmt += """
        AND window{cur_window} NOT IN (SELECT modelb FROM modeldist{prior_window}{cur_window}
        WHERE modela={{path[{prior_window}]}})""".format(cur_window=cur_window, prior_window=prior_window)

        if nwindows > 2:
            logging.debug("Next window query: %s", child_q_fmt)
            logging.debug("Edgescore index: %s", new_edgescore_index)

        with shared.ro_conn(start_db_file) as start_db:
            previous_path_list = list(start_db.cursor().execute(start_query))

        logging.debug("%s previous paths", len(previous_path_list))

        # Get single model scores
        score_q = "SELECT modelid, {model_score_columns} FROM model JOIN allmodel USING(modelid)".format(model_score_columns=", ".join(self.model_score_columns))
        with shared.ro_conn(score_db_file) as score_conn:
            score_gen = score_conn.cursor().execute(score_q)
            model_scores = {row[0]: dict(zip(self.model_score_columns, row[1:])) for row in score_gen}

        try:
            in_db_disk = apsw.Connection(in_db_file, flags=apsw.SQLITE_OPEN_READONLY)
        except apsw.CantOpenError:
            print in_db_file
            raise
        memcon = apsw.Connection(":memory:")
        with memcon.backup("main", in_db_disk, "main") as backup:
            backup.step()
        in_db_disk.close()
        in_db = memcon.cursor()

        mean = np.mean

        def gen_paths():
            for prev_path in previous_path_list:
                if nwindows == 2:
                    window_slice = 2
                else:
                    window_slice = cur_window
                # Only keep windows (remove score components)
                prev_path_windows = prev_path[:window_slice]
                row_dict = dict(zip(window_columns, prev_path_windows))
                # Prepare score components
                mean_prev_edgescore = prev_path[edgescore_index]
                di_list = [model_scores[modelid]['di'] for modelid in prev_path_windows]
                if nwindows == 2:
                    insert_dict = row_dict.copy()
                    insert_dict["edgescores"] = mean_prev_edgescore
                    insert_dict["nodescore"] = mean(di_list)
                    yield insert_dict
                else:
                    prev_edgescore = mean_prev_edgescore * n_prev_edges
                    current_path_q = child_q_fmt.format(path=prev_path_windows)
                    current_path_list = in_db.execute(current_path_q)
                    for new_path in current_path_list:
                        new_path_window = new_path[0]
                        new_path_edgescore = new_path[new_edgescore_index]
                        new_edgescore = (prev_edgescore + new_path_edgescore) / n_edges
                        new_di_list = list(di_list)
                        new_di_list.append(model_scores[new_path_window]['di'])
                        nodescore = mean(new_di_list)
                        insert_dict = row_dict.copy()
                        insert_dict[new_window] = new_path_window
                        insert_dict["edgescores"] = new_edgescore
                        insert_dict["nodescore"] = nodescore
                        yield insert_dict

        old_total_paths = total_paths
        try:
            out_db_conn = apsw.Connection(out_db_file, flags=apsw.SQLITE_OPEN_READWRITE | apsw.SQLITE_OPEN_CREATE)
        except apsw.CantOpenError:
            print out_db_file
            raise
        out_db_conn.cursor().execute(out_schema)
        # Outer transaction
        with out_db_conn:
            out_db_curs = out_db_conn.cursor()
            path_gen = gen_paths()
            for path_chunk in itertools.izip_longest(*[iter(path_gen)] * batchsize, fillvalue=None):
                # Inner transaction
                with out_db_conn:
                    for row_dict in path_chunk:
                        if row_dict is None:
                            continue
                        total_paths += 1
                        out_db_curs.execute(insert_sql, row_dict)

        # Cleanup
        memcon.close()
        out_db_conn.close()
        new_paths = old_total_paths - total_paths
        if not new_paths:
            raise FindPathsStepwiseError("No new paths were produced")

    @classmethod
    def commandline(cls, module_args=None):
        desc = """HELPDESCRIPTION"""
        a = argparse.ArgumentParser(description=desc)
        a.add_argument("complexid",
                       help="PDB identifier (e.g. 1ppeEI")
        a.add_argument("-d", "--directory", metavar="DIR",
                       help="Path to db files")
        a.add_argument("-n", "--nwindows", type=int, required=True,
                       help="Number of windows to calculate")
        a.add_argument("-b", "--batchsize", type=int, default=cls.default_batchsize,
                       help="Number of paths to check at a time")
        a.add_argument("-v", "--verbose", action="store_true",
                       help="Print debugging statements")

        args = a.parse_args(module_args)
        kwargs = vars(args)

        if kwargs.pop("verbose"):
            level = logging.DEBUG
        else:
            level = logging.INFO
        logging.basicConfig(level=level)

        c = cls(**kwargs)
        return c


if __name__ == "__main__":
    FindPathsStepwise.commandline()
    logging.info("Job complete after %s paths (%s skipped)", total_paths, skipped_paths)
