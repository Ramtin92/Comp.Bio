from __future__ import print_function

from Bio.KEGG import _write_kegg
from Bio.KEGG import _wrap_kegg


MODE_NORMAL = 0
MODE_RCLASS = 1

MODE_ORTHOLOGY = 3


parsing_mode = MODE_NORMAL

# Set up line wrapping rules (see Bio.KEGG._wrap_kegg)
name_wrap = [0, "",
             (" ", "$", 1, 1),
             ("-", "$", 1, 1)]
id_wrap = lambda indent: [indent, "", (" ", "", 1, 0)]
struct_wrap = lambda indent: [indent, "", ("  ", "", 1, 1)]


class Record(object):
    """Holds info from a KEGG Ligand/Compound record.

    Members:
    entry       The entry identifier.
    name        A list of the compund names.
    definition  The chemical definition for the compound
    equation
    rclass      A list of 2-tuples: (RC, RPAIR)
    pathway     A list of 3-tuples: (database, id, pathway)
    enzyme      A list of 2-tuples: (enzyme id, role)
    orthology      A list of 2-tuples: (knumber, kname)
    dblinks     A list of 2-tuples: (database, list of link ids)

    """

    def __init__(self):
        """__init___(self)

        Create a new Record.
        """
        self.entry = ""
        self.name = []
        self.definition = ""
        self.equation = ""
        self.rclass = []
        self.pathway = []
        self.enzyme = []
        self.orthology = []
        self.dblinks = []

    def __str__(self):
        """__str__(self)

        Returns a string representation of this Record.
        """
        return self._entry() + \
               self._name() + \
               self._definition() + \
               self._equation() + \
               self._rclass() + \
               self._pathway() + \
               self._enzyme() + \
               self._orthology() + \
               self._dblinks() + \
               "///"

    def _entry(self):
        return _write_kegg("ENTRY",
                           [self.entry])

    def _name(self):
        return _write_kegg("NAME",
                           [_wrap_kegg(l, wrap_rule=name_wrap)
                            for l in self.name])

    def _definition(self):
        return _write_kegg("DEFINITION",
                           [self.definition])

    def _equation(self):
        return _write_kegg("EQUATION",
                           [self.equation])

    def _rclass(self):
        s = []
        for entry in self.rclass:
            s.append(entry[0] + ": " + "  ".join(entry[1]) + "  ")
        return _write_kegg("RCLASS",
                           [_wrap_kegg(l, wrap_rule=id_wrap(0))
                            for l in s])

    def _pathway(self):
        s = []
        for entry in self.pathway:
            s.append(entry[0] + ": " + entry[1] + "  " + entry[2])
        return _write_kegg("PATHWAY",
                           [_wrap_kegg(l, wrap_rule=id_wrap(16))
                            for l in s])

    def _enzyme(self):
        s = ""
        for entry in self.enzyme:
            if entry[1]:
                t = entry[0] + " (" + entry[1] + ")"
            else:
                t = entry[0]
            s = s + t.ljust(16)
        return _write_kegg("ENZYME",
                            [_wrap_kegg(s, wrap_rule=id_wrap(0))])
    def _orthology(self):
        s = []
        for entry in self.orthology:
            s.append(entry[0] + ": " + "  ".join(entry[1]) + "  ")
        return _write_kegg("ORTHOLOGY",
                           [_wrap_kegg(l, wrap_rule=id_wrap(0))
                            for l in s])

    def _dblinks(self):
        s = []
        for entry in self.dblinks:
            s.append(entry[0] + ": " + " ".join(entry[1]))
        return _write_kegg("DBLINKS",
                           [_wrap_kegg(l, wrap_rule=id_wrap(9))
                            for l in s])


def parse(handle):
    """Parse a KEGG Ligan/Compound file, returning Record objects.

    This is an iterator function, typically used in a for loop.  For
    example, using one of the example KEGG files in the Biopython
    test suite,

    >>> with open("../../output/C00154.txt") as handle:
    ...     for record in parse(handle):
    ...         print("%s %s" % (record.entry, record.name[0]))
    ...
    C00023 Iron
    C00017 Protein
    C00099 beta-Alanine
    C00294 Inosine
    C00298 Trypsin
    C00348 Undecaprenyl phosphate
    C00349 2-Methyl-3-oxopropanoate
    C01386 NH2Mec

    """


    record = Record()
    for line in handle:
        if line[:3] == "///":
            yield record
            record = Record()
            continue
        if line[:12] != "            ":
            keyword = line[:12]
            if keyword != "RCLASS      " and keyword != "ORTHOLOGY   ":

                parsing_mode = MODE_NORMAL
        elif parsing_mode == MODE_RCLASS:
            keyword = "RCLASS      "
        elif parsing_mode == MODE_ORTHOLOGY:
            keyword = "ORTHOLOGY   "


        data = line[12:].strip()
        if keyword == "ENTRY       ":
            words = data.split()
            record.entry = words[0]
        elif keyword == "NAME        ":
            data = data.strip(";")
            record.name.append(data)
        elif keyword == "ENZYME      ":
            while data:
                column = data[:16]
                data = data[16:]
                if '(' in column:
                    entry = column.split()
                    enzyme = (entry[0], entry[1][1:-1])
                else:
                    enzyme = (column.strip(), "")
                record.enzyme.append(enzyme)
        elif keyword == "PATHWAY     ":
            if data[:5] == 'PATH:':
                path, map, name = data.split(None, 2)
                pathway = (path[:-1], map, name)
                record.pathway.append(pathway)
            # else:
            #     pathway = record.pathway[-1]
            #     path, map, name = pathway
            #     name = name + " " + data
            #     pathway = path, map, name
            #     record.pathway[-1] = pathway
        elif keyword == "DEFINITION  ":
            record.definition = data
        elif keyword == "EQUATION    ":
            record.equation = data
        elif keyword == "RCLASS      ":
            parsing_mode = MODE_RCLASS

            data_row = data.split()
            rcnumber = data_row[0]
            rpair = data_row[1]
            record.rclass.append((rcnumber, rpair))

        elif keyword == "ORTHOLOGY   ":
            parsing_mode = MODE_ORTHOLOGY

            data_row = data.split()
            knumber = data_row[0]
            kname = data_row[1]
            record.orthology.append((knumber, kname))

        elif keyword == "DBLINKS     ":
            if ":" in data:
                key, values = data.split(":")
                values = values.split()
                row = (key, values)
                record.dblinks.append(row)
            else:
                row = record.dblinks[-1]
                key, values = row
                values.extend(data.split())
                row = key, values
                record.dblinks[-1] = row


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
