#!/usr/bin/env python

import argparse
import logging
import os

import pandas as pd

logging.basicConfig(level=logging.DEBUG)


class SSError(RuntimeError):
    "Exception for SS problems"
    def __init__(self, message):
        super(SSError, self).__init__(message)

class ParseSs(object):

    fmt_str = "{index: >4} {aa} {ss}   {C:.3f}  {H:.3f}  {E:.3f}\n"
    pred_set = set(('C', 'E', 'H'))
    aa_set = set(('R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P',
                  'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'))

    @classmethod
    def ss_freq(cls, ss):
        return {k: 0.67 if ss == k else 0.15 for k in cls.pred_set}

    def __init__(self):
        """CONSTRUCTOR"""

    def read_file(self, method, infile, destdir=None, filebase=None):
        "File reading dispatch"
        basedir, filename = os.path.split(infile)
        if destdir is None:
            destdir = basedir
        self.basedir = destdir
        if filebase is None:
            filebase = filename.split(".")[0]
        self.filebase = filebase
        with open(infile, "r") as ih:
            filelines = list(ih)

        # Get method (raises AttributeError on invalid method)
        line_method = getattr(self, "read_%s" % method)

        lines = line_method(filelines)

        self.write(method, lines)

    def read_porter(self, filelines):
        "Read PORTER file (pasted from email to text)"
        lines = list()
        seq = ''
        pred = ''
        prev_line_seq = False
        for line in filelines:
            line = line.strip()
            # Comment lines contain lowercase
            if not line.isupper(): continue
            # Append to prediction if previous line was sequence and
            # all characters are SS types 
            if prev_line_seq and all(c in self.pred_set for c in line):
                pred += line
                prev_line_seq = False
            # Append to sequence if all characters are AA
            elif all(c in self.aa_set for c in line):
                seq += line
                prev_line_seq = True

        if len(seq) != len(pred):
            print len(seq), seq
            print len(pred), pred
            raise SSError("Length of sequence and prediction do not match")

        for x, (aa, ss) in enumerate(zip(seq, pred)):
            #print x + 1, s, p
            linedict = dict(index=x + 1,
                            aa=aa,
                            ss=ss)
            linedict.update(self.ss_freq(ss))
            lines.append(linedict)

        return lines

    def read_jpred(self, filelines):
        "Read Jpred .concise output file"
        match_key_dict = dict(jnetpred='ss',
                              #QUERY='aa',
                              #align1='aa',
                              JNETPROPE='E',
                              JNETPROPH='H',
                              JNETPROPC='C')
        # Avoid matching 'align19'
        match_key_dict['align1;'] = 'aa'
        data_dict = dict()
        for line in filelines:
            if not line.strip(): continue
            # Remove stuff after trailing comma
            parts = [p.strip() for p in line.split(",") if p.strip()]
            label, zero = parts[0].split(":")
            key = [v for k, v in match_key_dict.iteritems() if k in label]
            if not key: continue
            if len(key) > 1:
                raise SSError("Ambiguous label %s" % label)
            key = key[0]
            parts[0] = zero
            data_dict[key] = parts

        if len(set([len(v) for v in data_dict.values()])) != 1:
            print len(data_dict['seq']), "".join(data_dict['seq'])
            print len(data_dict['pred']), "".join(data_dict['pred'])
            raise SSError("Length of sequence and prediction do not match")

        data_df = pd.DataFrame(data_dict)
        for col in "C", "H", "E":
            data_df[col] = pd.to_numeric(data_df[col])
        data_df = data_df.dropna()
        data_df.loc[:, 'ss'] = data_df['ss'].replace('-', 'C')
        data_df.loc[:, 'index'] = range(1, len(data_df) + 1)

        return data_df.to_dict("records")

    def read_generic(self, filelines):
        "Fragile: assumes first line is sequence and second is prediction"
        lines = list()
        counter = 0
        for line in filelines:
            line = line.strip()
            if not line or line[0] == "#": continue
            counter += 1
            text = line.split()[-1]
            if counter == 1:
                seq = text
            elif counter == 2:
                pred = text

        if len(seq) != len(pred):
            print len(seq), "".join(seq)
            print len(pred), "".join(pred)
            raise SSError("Length of sequence and prediction do not match")

        for x, (aa, ss) in enumerate(zip(seq, pred)):
            if ss == "-":
                ss = "C"
            linedict = dict(index=x + 1,
                            aa=aa,
                            ss=ss)
            linedict.update(self.ss_freq(ss))
            lines.append(linedict)

        return lines

    def read_sspro(self, filelines):
        lines = list()
        seq_index = None
        pred_index = None
        for x, line in enumerate(filelines):
            if line.startswith("Amino Acids:"):
                seq_index = x + 1
            elif line.startswith("Predicted Secondary Structure"):
                pred_index = x + 1
            elif x == seq_index:
                seq = line.strip()
            elif x == pred_index:
                pred = line.strip()

        if len(seq) != len(pred):
            print len(seq), "".join(seq)
            print len(pred), "".join(pred)
            raise SSError("Length of sequence and prediction do not match")

        for x, (aa, ss) in enumerate(zip(seq, pred)):
            linedict = dict(index=x + 1,
                            aa=aa,
                            ss=ss)
            linedict.update(self.ss_freq(ss))
            lines.append(linedict)

        return lines

    def write(self, method, lines):
        "Write currently loaded file and clear"
        with open(os.path.join(self.basedir, "%s.%s.ss2" % (self.filebase, method)), "w") as oh:
            oh.write("# PSIPRED VFORMAT (PSIPRED V3.3)\n\n")
            for line_dict in lines:
                oh.write(self.fmt_str.format(**line_dict))

    @classmethod
    def commandline(cls, module_args=None):
        desc = """Parser for secondary structure prediction formats"""
        a = argparse.ArgumentParser(description=desc)
        a.add_argument("--generic",
                       help="Generic format (first line: amino acid sequence, second line: 3-class SS pred)")
        a.add_argument("--porter",
                       help="PORTER email output")
        a.add_argument("--jpred",
                       help="Jpred .concise output")
        a.add_argument("--sspro",
                       help="SSpro email output")

        args = a.parse_args(module_args)
        c = cls()

        if args.porter:
            c.read_file("porter", args.porter)
        if args.jpred:
            c.read_file("jpred", args.jpred)
        if args.generic:
            c.read_file("generic", args.generic)
        if args.sspro:
            c.read_file("sspro", args.sspro)

        return c


if __name__ == "__main__":
    ParseSs.commandline()
