"""
Collection of old parsers.
"""

__license__ = """
    Copyright 2008-2012,2016 Stephane De Mita, Mathieu Siol

    This file is part of EggLib.

    EggLib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    EggLib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with EggLib.  If not, see <http://www.gnu.org/licenses/>.
"""

from .. import _eggwrapper, _interface
import string, StringIO, re

########################################################################

def from_clustal(string):

    """
    Imports a clustal-formatted alignment. The input format is the one
    generated and used by CLUSTALW (see `<http://web.mit.edu/meme_v4.9.0/doc/clustalw-format.html>`_).

    :param string: input clustal-formatted sequence alignment.
    :return: A new :class:`.Align` instance.

    .. versionchanged:: 3.0.0
        Renamed (previous name was :func:`.aln2fas`). Input argument is
        a string rather than a file.
    """

    stream = StringIO.StringIO(string)

    # the first line must start by CLUSTAL W or CLUSTALW
    line = stream.readline()
    c = 1
    if (line[:8] != 'CLUSTALW' and
        line[:9] != 'CLUSTAL W' and
        line[:8] != 'CLUSTAL '): raise ValueError, 'invalid clustalw-formatted string (line: {0})'.format(c)

    # start reading blocks
    cnt = _eggwrapper.DataHolder(False)
    line = stream.readline()

    # initialize names list
    names = []

    while True:
        
        # skip empty lines
        while line != '' and line.strip() == '':
            line = stream.readline()
            c += 1

        # detect end of file
        if line == '': break

        # read a block
        block_length = None

        # read sequences
        while True:

            # conservation line
            if line[0] == ' ':
                if block_length == None: raise ValueError, 'invalid clustalw-formatted string (line: {0})'.format(c)
                line = stream.readline()
                c += 1
                if line == '': break # end of file
                if line.strip() != '': raise ValueError, 'invalid clustalw-formatted string (line: {0})'.format(c)
                break

            # reads sequence
            bits = line.split()

            # sequence line
            if len(bits) > 3:
                raise ValueError, 'invalid clustalw-formatted string (line: {0})'.format(c)

            name = bits[0]
            seq = bits[1]
            if block_length == None:
                block_length = len(seq)
            elif block_length != len(seq):
                raise ValueError, 'invalid clustalw-formatted string (line: {0})'.format(c)

            # adds the sequence to the container (new sequence)
            if name not in names:
                pos = len(names)
                cnt.set_nsam_i(pos + 1)
                cnt.set_name_i(pos, name)
                cnt.set_nsit_i(pos, len(seq))
                for i, v in enumerate(seq): cnt.set_i(pos, i, ord(v))
                names.append(name)

            # adds the sequence (continuing old sequence)
            else:
                pos = names.index(name)
                cur = cnt.get_nsit_i(pos)
                cnt.set_nsit_i(pos, cur + len(seq))
                for i, v in enumerate(seq): cnt.set_i(pos, cur+i, ord(v))

            if len(bits) == 3:
                try:
                    i = int(bits[2])
                except ValueError:
                    raise ValueError, 'invalid clustalw-formatted string (line: {0})'.format(c)
                if cnt.get_nsit_i(pos) != i:
                    raise ValueError, 'invalid clustalw-formatted string (line: {0})'.format(c)

            # checks next line
            line = stream.readline()
            c += 1
            if line == '':
                raise ValueError, 'invalid clustalw-formatted string (line: {0})'.format(c)
                
            # empty conservation line is caught by this line
            if line.strip() == '': break

    if not cnt.is_equal(): raise ValueError, 'invalid clustalw-formatted string (unequal sequences)'
    return _interface.Align._create_from_data_holder(cnt)

########################################################################

_staden_table = string.maketrans(' -*', '?N-')

########                                                        ########

def from_staden(string, delete_consensus=True):

    """
    Import the output file of the GAP4 program of the
    `Staden <http://staden.sourceforge.net/>`_ package.

    The input file should have been generated from a contig alignment by
    the GAP4 contig editor, using the command "dump contig to file". The
    sequence named ``CONSENSUS``, if present, is automatically removed
    unless the option *delete_consensus* is ``False``.

    Staden's default convention is followed:

    * ``-`` codes for an unknown base and is replaced by ``N``.
    * ``*`` codes for an alignment gap and is replaced by ``-``.
    * ``.`` represents the same sequence than the consensus at that
      position.
    * White space represents missing data and is replaced by ``?``.
    
    .. versionadded:: 2.0.1
        Add argument *delete_consensus*.
        
    .. versionchanged:: 2.1.0
        Read from string or fname.

    .. versionchanged:: 3.0.0
        Renamed :func:`.from_staden`. Only string input is supported
        now.
    """

    # get shift from the first CONSENSUS line
    mo = re.search('(        CONSENSUS +)[A-Za-z\-\*]', string)
    if mo == None: raise ValueError, 'invalid staden contig dump file'
    shift = len(mo.group(1))

    # split lines and identify blocks (based on empty lines)
    lines = string.splitlines()
    empty_lines = [i for i, v in enumerate(lines) if v == '']
    empty_lines.insert(0, -1) # emulate a white line immediately before first line
    empty_lines.append(len(lines))
    blocks = [lines[empty_lines[i]+2 : empty_lines[i+1]] for i in xrange(len(empty_lines) - 1)]
        # +1 to read after each blank line
        # +2 to skip the first line of each block

    # initialize variables
    align = _eggwrapper.DataHolder(False)
    ns = 0
    currpos = 0
    IDs = []

    # process all blocks
    for block in blocks:
        for line in block:
            ID = line[:7].strip()
            name = line[8:shift].strip()
            sequence = line[shift:]
            sequence = sequence.translate(_staden_table)
            block_length = len(sequence)

            if ID not in IDs:
                ns += 1
                align.set_nsam_i(ns)
                align.set_name_i(ns-1, name)
                align.set_nsit_i(ns-1, currpos)
                for i in xrange(currpos): align.set_i(ns-1, i, 63)
                pos = ns - 1
                IDs.append(ID)
            else:
                pos = IDs.index(ID)

            align.set_nsit_i(pos, currpos + len(sequence))
            for i,v in enumerate(sequence): align.set_i(pos, currpos + i, ord(v))

        currpos += block_length

    # equalize
    for i in xrange(ns):
        n = align.get_nsit_i(i)
        align.set_nsit_i(i, currpos)
        for j in xrange(n, currpos): align.set_i(i, j, 63)
    align.set_is_matrix(True)

    # undot
    idx = IDs.index('')
    for i in xrange(ns):
        if i != idx:
            for j in xrange(currpos):
                if align.get_i(i, j) == 46:
                    align.set_i(i, j, align.get_i(idx, j))

    # remove consensus
    if delete_consensus:
        align.del_sample_i(idx)

    # return
    return _interface.Align._create_from_data_holder(align)

########################################################################

def from_genalys(string):
    
    """
    Converts Genalys-formatted sequence alignment files to fasta. This
    function imports files generated through the option *Save SNPs* of
    Genalys 2.8.

    :param string: input data as a Genalys-formatted string.
    :return: An :class:`.Align` instance.

    .. versionchanged:: 3.0.0
        Renamed :func:`.from_genalys`. Only string input is supported
        now.
    """

    stream = StringIO.StringIO(string)

    insertions = []
    flag = False

    for line in stream:
        line = line.split("\t")

        if len(line) > 1 and line[0] == "Polymorphism":
            flag = True

        if len(line) > 1 and line[0] == "IN" and flag:
            insertions.extend(line[1].split("/"))

    if len(insertions) > 0:
        tp = insertions[0].split("_")
        if len(tp) == 1:
            tp = tp[0].split(".")
            if len(tp) == 1:
                tp.append("1")
        finsertions = [tp]

    for i in insertions:
        i = i.split("_")
        if len(i) == 1:
            tp = tp[0].split(".")
            if len(tp) == 1:
                i.append("1")
        if i[0] != finsertions[-1][0]:
            finsertions.append(i)
        finsertions[-1][1] = i[1]
    
    if len(insertions) > 0:
        insertions = finsertions

    stream.close()
    stream = StringIO.StringIO(string)

    names = []
    sequences = []
    maxlen = 0

    for line in stream:
        line = line.split("\t")
            
        if len(line) > 1:
            bidon = re.match(".+\.ab1$", line[1])
            if bidon != None:
                names.append(line[1])
                sequences.append("")
                index = 6
                for i in xrange(10):
                    if line[i] == "F" or line[i] == "R":
                        index = i+1
                        break
                if line[index] != "":
                    beginning = int(line[index]) - 1
                    for i in insertions:
                        if int(i[0]) <= beginning:
                            beginning = beginning + int(i[1])
                        else:
                            break
                    for i in xrange(beginning):
                        sequences[-1]= sequences[-1] + "?"
                sequences[-1] = sequences[-1] + line[-1].rstrip("\n")
                if len(sequences[-1]) > maxlen:
                    maxlen = len(sequences[-1])

    data = _eggwrapper.DataHolder(True)
    data.set_nsam_i(len(sequences))
    data.set_nsit(maxlen)

    for i in xrange(len(sequences)):
        data.set_name_i(i, names[i])
        sequences[i] = sequences[i].replace("_", "-")
        for j, v in enumerate(sequences[i]): data.set_i(i, j, ord(v))
        for j in xrange(len(sequences[i]), maxlen): data.set_i(i, j, 63)

    return _interface.Align._create_from_data_holder(data)

########################################################################

def get_fgenesh(string, locus='locus'):
    
    """
    Imports fgenesh output.

    :param fname: a string containing fgenesh ouput.
    :parma locus: locus name.
    :return: A :class:`list` of ``gene`` and ``CDS`` features
        represented by dictionaries. Note that 5' partial features
        might not be in the appropriate frame and that it can be
        necessary to add a ``codon_start`` qualifier.

    .. versionchanged:: 3.0.0
        Input as string. Added *locus* argument.
    """

    # supports for mac/windows files
    string = '\n'.join(string.splitlines())

    # gets the feature table
    try:
        data_sub = string.split('   G Str   Feature   Start        End    Score           ORF           Len\n')[1].split('Predicted protein(s):\n')[0]
    except IndexError:
        raise ValueError, 'invalid fgenesh format'
    data_sub = data_sub.split('\n\n')

    # edit
    del data_sub[-1]
    data_sub[0]= '\n'.join(data_sub[0].split('\n')[1:])

    # iteratively grabs the features
    features = {}
    for i in data_sub:
        pos = []
        start = 1
        rank = '---'
        strand = '---'
        for j in i.split('\n'):
            a = re.search(' ?[0-9]+ ([+|-])      (TSS|PolA) +([0-9]+)', j)
            b = re.search(' ?([0-9]+) ([+|-]) + ([0-9])+ CDS(o|f|i|l) +([0-9]+) - +([0-9]+) +[-\.0-9]+ + ([0-9]+)', j)
            if b:
                if b.group(3) == "1":
                    if int(b.group(5)) == int(b.group(7)): start= 1
                    elif int(b.group(5)) == (int(b.group(7))-1): start= 2
                    elif int(b.group(5)) == (int(b.group(7))-2): start= 3
                    else: raise ValueError, 'invalid fgenesh format'
                pos.append( [int(b.group(5))-1, int(b.group(6))-1 ] )
                rank = b.group(1)
                if b.group(2) == '+': strand = 'plus'
                else: strand = 'minus'

        features['cds'+rank] ={
            'gene': locus+'_'+rank,
            'strand': strand,
            'pos': pos,
            'type': 'CDS',
            'note': 'fgenesh prediction'
        }

        features['gene'+rank] ={
            'gene': locus+'_'+rank,
            'strand': strand,
            'pos': [[ pos[0][0], pos[-1][1] ]],
            'type': 'gene',
            'note': 'fgenesh prediction'
        }
        
    # gets the sequence section
    try:
        data_sub = string.split('   G Str   Feature   Start        End    Score           ORF           Len\n')[1].split('Predicted protein(s):\n')[1].split('>')
    except IndexError:
        raise ValueError, 'invalid fgenesh format'
    del data_sub[0]

    if ( (2*len(data_sub) != len(features)) and
           (len(data_sub) != len(features)) ) : raise ValueError, 'invalid fgenesh format'

    # returns the sequences as a table
    return list(features.itervalues())
