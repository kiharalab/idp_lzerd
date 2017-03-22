# Copyright 2013-2016 Lenna X. Peterson. All rights reserved.

from contextlib import contextmanager
import errno
import inspect
import itertools
import math
import os
import re
import StringIO

import apsw
import pandas as pd

script_dir = os.path.abspath(os.path.dirname(inspect.getfile(inspect.currentframe())))
# Root directory of package
ROOTDIR = os.path.normpath(os.path.join(script_dir, os.pardir))


class IDPError(RuntimeError):
    """Exception for IDP-LZerD"""
    def __init__(self, message):
        super(IDPError, self).__init__(message)

#### SQL FUNCTIONS ####
@contextmanager
def ro_conn(dbfile):
    try:
        with apsw.Connection(dbfile, flags=apsw.SQLITE_OPEN_READONLY) as conn:
            yield conn
    except apsw.CantOpenError:
        print dbfile
        raise

@contextmanager
def write_conn(dbfile):
    try:
        with apsw.Connection(dbfile, flags=apsw.SQLITE_OPEN_READWRITE) as conn:
            yield conn
    except apsw.CantOpenError:
        print dbfile
        raise

@contextmanager
def new_conn(dbfile):
    try:
        with apsw.Connection(dbfile, flags=apsw.SQLITE_OPEN_CREATE | apsw.SQLITE_OPEN_READWRITE) as conn:
            yield conn
    except apsw.CantOpenError:
        print dbfile
        raise

def db_to_pandas(select, dbf, **kwargs):
    """
    Open database and load select into pandas DataFrame.

    :param select: SQL select query
    :param dbf: database file
    :param kwargs: additional arguments to pandas.read_sql_query

    :returns: pandas.DataFrame
    """
    with ro_conn(dbf) as conn:
        return conn_to_pandas(select, conn, **kwargs)

def conn_to_pandas(select, conn, **kwargs):
    """
    Load the results of a select into a pandas DataFrame.

    :param select: SQL select query
    :type select: str
    :param conn: database connection
    :type conn: apsw.Connection
    :param kwargs: additional arguments to pandas.read_sql_query

    :returns: pandas.DataFrame
    """
    try:
        return pd.read_sql_query(select, conn, **kwargs)
    except apsw.ExecutionCompleteError:
        return pd.DataFrame()

def create_insert_statement(tablename, columns):
    """
    Create insert statement.
    e.g. func("mytable", ["foo", "bar"])
    INSERT INTO mytable (foo, bar) VALUES (:foo, :bar)

    :param tablename: name of table to update
    :type tablename: str
    :param columns: columns to set
    :type columns: list
    """
    insert_str = "INSERT INTO {tablename} ({columns}) VALUES ({bindings})"
    return insert_str.format(tablename=tablename,
                             columns=", ".join(columns),
                             bindings=", ".join([":" + v for v in columns]))

def create_update_statement(tablename, columns, where):
    """
    Create update statement.
    e.g. func("mytable", ["foo", "bar"], ["id"])
    UPDATE mytable SET foo=:foo, bar=:bar WHERE id=:id

    :param tablename: name of table to update
    :type tablename: str
    :param columns: columns to set
    :type columns: list
    :param where: where columns to match
    :type where: list

    :returns: str
    """
    if isinstance(where, basestring):
        where = [where]
    if isinstance(columns, basestring):
        columns = [columns]
    update_fmt = "UPDATE {tablename} SET {values} WHERE {where}"
    update_stmt = update_fmt.format(tablename=tablename,
                                    values=", ".join(["{0}=:{0}".format(c) for c in columns]),
                                    where=" AND ".join(["{0}=:{0}".format(w) for w in where]))
    return update_stmt

#### FILE SYSTEM FUNCTIONS ####
class CHDIR(object):
    def __init__(self, dirname):
        self.dirname = dirname

    def __enter__(self):
        self.old_dir = os.getcwd()
        os.chdir(self.dirname)

    def __exit__(self, exception_type, exception_value, traceback):
        os.chdir(self.old_dir)


def missing(*files):
    "Check whether files are missing or empty"
    return any(not os.path.isfile(f) or not os.path.getsize(f)
               for f in files)

def mkdir_p(path):
    """
    Pythonic replacement for shell `mkdir -p`
    Creates multiple directories and ignores error caused by existing
    """
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def silent_remove(path):
    """
    Remove file, ignore if missing.
    """
    try:
        os.remove(path)
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise

#### PDB FUNCTIONS ####
_h_re = re.compile(r"[123 ]*H.*")

def strip_h(filename):
    """
    Strip hydrogens from PDBfile
    NB loads entire file into memory, use at your own risk
    """
    with open(filename, "r") as ih:
        ret = StringIO.StringIO("\n".join(r for r in ih
                                          if r[:4] != "ATOM"
                                          or not _h_re.match(r[12:16])))
    return ret

three_to_one = {
    'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'GLU': 'E', 'GLN': 'Q',
    'ASP': 'D', 'ASN': 'N', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F', 'TYR': 'Y',
    'ARG': 'R', 'LYS': 'K', 'SER': 'S', 'THR': 'T', 'MET': 'M', 'ALA': 'A',
    'GLY': 'G', 'PRO': 'P', 'CYS': 'C'}

one_to_three = {o: t for t, o in three_to_one.iteritems()}

#### SCORING FUNCTIONS ####
index_names = dict(model="modelindex",
                   fragment="fragmentindex",
                   relax="pathsid",)

index_re = dict(model=re.compile(r"model(\d+)\D*"),
                fragment=re.compile(r'frag_(\d\d\d)_model_h.pdb'),
                relax=re.compile(r"\D*(\d+)\D*"))

def df_extract_index(df, kind, colname="model"):
    indexname = index_names[kind]
    regex = index_re[kind]
    df[indexname] = df.pop(colname).str.extract(regex, expand=False)
    df[indexname] = df[indexname].astype(int)
    return df

itscore_file = "scores.itscore"
itscore_header = ("model", "itscore")
itscore_types = (str, float)

def read_itscore(filename, kind=None):
    kwargs = dict(dtype=zip(itscore_header, itscore_types),
                  delim_whitespace=True,
                  names=itscore_header)
    mydf = pd.read_table(filename, **kwargs)
    if kind:
        mydf = df_extract_index(mydf, kind=kind)
    return mydf

goap_file = "goap_score.txt"
goap_header = ("index", "model", "goap", "dfire", "extra")
goap_types = (str, str, float, float, float)

def read_dfire_from_goap(filename, all=False, kind=None):
    kwargs = dict(dtype=zip(goap_header, goap_types),
                  names=goap_header,
                  delim_whitespace=True)
    if not all:
        kwargs['usecols'] = ['model', 'dfire']
    mydf = pd.read_csv(filename, **kwargs)
    if kind:
        mydf = df_extract_index(mydf, kind=kind)
    return mydf

#### IDP-LZERD FUNCTIONS ####
def create_windows(length):
    """Compute and yield fragments"""
    last_loop = False
    for start, end in itertools.izip(itertools.count(1, 6), itertools.count(9, 6)):
        position = start
        if end >= length:
            last_loop = True
            end = length
            position = end - 8
            start = int(6.0 * math.floor((position + 4.4) / 6.0) + 1.0)
            if position < 1:
                raise IDPError("Length %s is too short" % length)
        yield dict(position=position, res_start=start, res_end=end, windowindex=start)
        if last_loop: break

def add_model_path(windowrows, model_db_file):
    """
    Add model path to windows
    """
    model_key = "modelid"
    window_key = "window"
    path_key = "pathsid"
    model_query = """SELECT modelid, modelindex, fragmentindex, window_wd
    FROM model JOIN allmodel USING(modelid)
    JOIN window USING(windowindex)
    """
    model_data = db_to_pandas(model_query, model_db_file, index_col=model_key)
    model_data['path'] = model_data.apply(lambda x: os.path.join(x['window_wd'], str(x['fragmentindex']), "decoys", "model%s.pdb" % x['modelindex']), 1)

    # Associate paths with modelids
    value_vars = [x for x in windowrows.columns.values.tolist() if x.startswith(window_key)]
    melted = pd.melt(windowrows, id_vars=path_key, value_vars=value_vars, var_name=window_key, value_name=model_key)
    # Double bracketed column selection keeps result as DataFrame
    melted = melted.merge(model_data.loc[:, ['path']], how="left", left_on=model_key, right_index=True)
    def fake_agg(g):
        T = tuple(g)
        if len(T) > 1:
            return T
        else:
            return T[0]
    newwindowrows = pd.pivot_table(melted.loc[:, [path_key, window_key, 'path']], index=path_key, columns=window_key, aggfunc=fake_agg)
    # levels are ('path', 'windowN')
    newwindowrows.columns = newwindowrows.columns.get_level_values(1)
    # will put index 'pathsid' into column
    newwindowrows = newwindowrows.reset_index()
    return newwindowrows


def load_config():
    # Parse paths
    path_file = os.path.join(ROOTDIR, "PATHS.ini")
    if not os.path.isfile(path_file):
        raise IDPError("Could not open PATHS.ini")
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
                    print "Could not interpret config line %s" % line

    required_keys = ('lzerd_path', 'rosetta_path', 'nr_path', 'blastpgp_exe')

    #if 'lzerd_path' not in config:
        #raise IDPError("'{0}' did not contain required key 'lzerd_path'".format(path_file))
    errors = 0
    for key in required_keys:
        if key not in config:
            errors += 1
            raise IDPError("'{0}' did not contain required key '{1}'".format(path_file, key))
    if not errors:
        return config
