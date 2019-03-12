"""
The :func:`.blast` function.
"""

__license__ = """
    Copyright 2009-2017 Stephane De Mita, Mathieu Siol, Thomas Coudoux

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

import subprocess, re, os
import xml.dom.minidom as MM
import xml.etree.ElementTree as ET
from .. import _interface
from ..io import _fasta
import _utils

class Hsp_data(object):  #High-scoring Segment Pair
    def __init__(self):
        self._num = None
        self._id = None
        self._def = None
        self._accession = None
        #self._len = None
        self._bit_score = None
        self._score = None
        self._evalue = None
        self._query_start = None
        self._query_end = None
        self._hit_start = None
        self._hit_end = None
        self._query_frame = None
        self._hit_frame = None
        self._identity = None
        self._positive = None
        self._gaps = None
        self._align_len = None
        self._qseq = None
        self._hseq = None
        self._midline = None

    num = property(lambda self: self._num, doc='index of the hit sequence')
    Id = property(lambda self: self._id, doc='ID of the hit sequence')
    Def =  property(lambda self: self._def, doc='definition line of subject')
    accession = property(lambda self: self._accession, doc='the accession of the matched database sequence')
    #Len = property(lambda self: self._len, doc='')
    bit_score = property(lambda self: self._bit_score, doc='The higher the bit-score, the better the sequence similarity')
    score = property(lambda self: self._score, doc='numerical value that describes the overall quality of an alignment.')
    evalue = property(lambda self: self._evalue, doc='E-value is the number of expected hits of similar quality (score) that could be found just by chance.')
    query_start = property(lambda self: self._query_start, doc='start of alignment in query')
    query_end = property(lambda self: self._query_end, doc='end of alignment in query')
    hit_start = property(lambda self: self._hit_start, doc='start of alignment in subject')
    hit_end = property(lambda self: self._hit_end, doc='end of alignment in subject')
    query_frame = property(lambda self: self._query_frame, doc='Query frame')
    hit_frame = property(lambda self: self._hit_frame, doc='Subject frame')
    identity = property(lambda self: self._identity, doc='the highest percent identity')
    positive = property(lambda self: self._positive, doc='Percentage of positive-scoring matches')
    gaps = property(lambda self: self._gaps, doc='Total number of gaps')
    align_len = property(lambda self: self._align_len, doc='alignment length')
    qseq = property(lambda self: self._qseq, doc='Aligned part of query sequence')
    hseq = property(lambda self: self._hseq, doc='Aligned part of subject sequence')
    midline = property(lambda self: self._midline, doc='formating middle line of the sequence aligment')


class Alignment(object):
    """ 

    """

    def __init__(self):
        self._iter_num = None
        #self._query_ID = None
        self._query_def = None
        self._query_len = None
        self._hit = []
        self._statistics = {}

    def _append_hit(self, value):
        self._hit.append(value)

    def __getitem__(self, index):
        if index > len(self._hit) or index < 0:
           raise ValueError, 'invalid index'
        else:
           return self._hit[index]

    def __iter__(self):
        for i in xrange (len(self._hit)):
            yield self._hit[i]

    def __len__(self):
        return len(self._hit)

    iter_num = property(lambda self: self._iter_num, doc='iteration number')
    query_ID = property(lambda self: self._query_ID, doc='definition line of subject')
    query_def = property(lambda self: self._query_def, doc='Definition line of query')
    query_len = property(lambda self: self._query_len, doc='length of query sequence')


class Blast(object):
    """ 
    Blast holds all Blast's informations, from an xml output.
    This class contains three variables wich representing all parts of a Blast output.
    :Db: name of the database used to make the Blast alignement (Contains subject sequence).
    :Parameters: Dictionnary containing all information about the configuration of the current Blast
    :Alignments: list containing all informations about blast replace gap by another characterone or multiple aligments which occurred. One
                 item in this list represent an aligment between a query sequence and differents subject
                 sequences. This list, contains object of the class Alignment.

    :class:`~._wrappers.Blast` instances are iterable by the special 
    :meth:`~._wrappers.Blast.__iter__` (with support 'for' statement).
    The iteration is done on alignments in the blast object. Each 
    itteration return a :class:`~._wrappers.Alignment` instance.
    :class:`~._wrappers.Blast` instances are slices objects (usable by the next
    command: self[index]). The method :meth:`~._wrapper.Blast.__getitem__` allows 
    the access by index. This method gets the alignement at given index. Negatives
    values are not allowed, else returns a 'ValueError'.

    This class is hiding in the module blast, so cannot be created without call one of the function blast
    """

    def __init__(self): 
        self._db = None
        self._parameters = {}
        self._alignements = []

    def __len__(self):
        """
        get the number of sites in the current :class:`~._vcf.BedWindow`
        instance
        this method is usable in this way: len(self)
        """
        return len(self._alignements)

    def __getitem__(self, index):
        """
        Get the site at given index (used in expression: self[indice]). 
        """ 
        if index > len(self._alignements) or index < 0:
           raise ValueError, 'invalid index'
        else:
           return self._alignements[index]

    def __iter__(self):
        """
        Implementation of the iterable interface of the :class:`~._wrappers.Blast`
	The iteration is done on alignments. 
        So the :class:`~._wrapper.Blast` instances are iterable and can
	be use with support 'for' 
        """
        for i in xrange (len(self._alignements)):
            yield self._alignements[i]


    def _append_align(self, value):
        """
        This method allows to append the list aligment of 
        """
        if not isinstance(value, Alignment): raise TypeError, 'The value {0} is not an object of the class alignement'
        self._alignements.append(value)

    db = property(lambda self: self._db, doc='Name of the database')
    parameters = property(lambda self: self._parameters, doc='information about the BLAST parameters.')


def _xml_config(fname, output='tmp.xml'):
    """ 
    This method allow to format a xml files according the tag. To avoid the indentation problems
    :fname: xml file
    :output: output, xml file formated
    """
    if not os.path.isfile(fname): raise  ValueError, 'The file {0} doesn`t exist'.format(source)
    string = open(fname, 'r').read()
    string = re.sub('\<\?.*\n|\<\!.*\n', '', string )
    xml_str = string.replace('\n', '')
    re.sub('\>( *?)\<', '', xml_str)
    xml_conf = MM.parseString(xml_str).toprettyxml(indent=" ")
    with open(output, "w") as f:
        f.write(xml_conf)


def _value_format(value):
    """
    Convert a digit string in the correct type: Int, Float, negative of positive.
    :value: a string containing digit value.
    """

    if value.replace('-','').isdigit():
	if '-' in value: return -int(value[1:])
        else: return int(value);
    elif '.' in value and value.replace('.','').replace('e','').replace('-','').isdigit():
	if '-' in value: return -float(value[1:])
        else: return float(value);
    else:
        return value;
 

def _get_aligmements(node):
    """ 
    This method allows to holds information about hits in the alignments section and 
    information about hsps in an alignment hit. Each database search hit comprises 1
    or more HSPs. Each Hsp has a corresponding score, significance level and the 
    fragment count, N...
     
    """
    if not isinstance(node, ET.Element): raise ValueError, 'invalid value for argument {0}'.format(node)
    Hits = Alignment()
    Hits._iter_num = _value_format(node.find('Iteration_iter-num').text)
    Hits._query_ID = node.find('Iteration_query-ID').text
    Hits._query_def = node.find('Iteration_query-def').text
    Hits._query_len = _value_format(node.find('Iteration_query-len').text)
    if len(node.find('Iteration_hits')) > 0:
        for Hit in node.find('Iteration_hits'):
            hsp = Hsp_data()
            hsp._num = _value_format(Hit.find('Hit_num').text)
            hsp._id = Hit.find('Hit_id').text
            hsp._def =  Hit.find('Hit_def').text
            hsp._accession = _value_format(Hit.find('Hit_accession').text)
            
            Hit_hsps = Hit.find('Hit_hsps')
            Hsp = Hit_hsps.find('Hsp')
            
            hsp._len =  _value_format(Hsp.find('Hsp_align-len').text)
            hsp._bit_score = _value_format(Hsp.find('Hsp_bit-score').text)
            hsp._score = _value_format(Hsp.find('Hsp_score').text)
            hsp._evalue = _value_format(Hsp.find('Hsp_evalue').text)
            hsp._query_start = _value_format(Hsp.find('Hsp_query-from').text)
            hsp._query_end = _value_format(Hsp.find('Hsp_query-to').text)
            hsp._hit_start = _value_format(Hsp.find('Hsp_hit-from').text)
            hsp._hit_end = _value_format(Hsp.find('Hsp_hit-to').text)
            hsp._query_frame = _value_format(Hsp.find('Hsp_query-frame').text)
            hsp._hit_frame =  _value_format(Hsp.find('Hsp_hit-frame').text)
            hsp._identity =  _value_format(Hsp.find('Hsp_identity').text)
            hsp._positive =  _value_format(Hsp.find('Hsp_positive').text)
            hsp._gaps =  _value_format(Hsp.find('Hsp_gaps').text)
            hsp._align_len = _value_format(Hsp.find('Hsp_align-len').text)
            hsp._qseq = Hsp.find('Hsp_qseq').text
            hsp._hseq = Hsp.find('Hsp_hseq').text
            hsp._midline = Hsp.find('Hsp_midline').text
            

        Itr_stats = node.find('Iteration_stat')
        Statistics = Itr_stats.find('Statistics')
        for stat in Statistics:
            tag = stat.tag.replace('Statistics_','')
            value = _value_format(stat.text)
            Hits._statistics[tag] = value
 
        Hits._append_hit(hsp)

    return Hits

def get_data_xml(fname):
    """ 
    This function allows to gets information of an XMl blast ouput and saves them in a instance of the :class:`~._wrappers.Blast`
    """
    _xml_config(fname)
    BLAST = ET.parse('tmp.xml')    
    root = BLAST.getroot()
    if not isinstance(root, ET.Element): raise ValueError, 'invalid value for argument {0}'.format(root)
    data = root.getchildren()

    BOp_params = root.find('BlastOutput_param')
    Parameters = BOp_params.find('Parameters')
    BR = Blast()
    BR._db = root.find('BlastOutput_db').text

    for param in Parameters:
        tag = param.tag.replace('Parameters_','')
        value = _value_format(param.text)
        BR._parameters[tag] = value

    BOp_itrs = root.find('BlastOutput_iterations')
    tmp_al=()
    for itr in BOp_itrs:    
        alignment = _get_aligmements(itr)
        BR._append_align(alignment)
    return BR


#Contains all commons arguments between the Blast tools
_common_args = { 
    'best_hit_overhang':        (float, lambda x: x>=0.0 and x<=0.5, lambda x: True if set(x).isdisjoint(['-culling_limit']) else False),
    'best_hit_score_edge':        (float, lambda x: x>=0.0 and x<=0.5, lambda x: True if set(x).isdisjoint(['-culling_limit']) else False),
    'culling_limit':        (int, lambda x: x>=0, lambda x: True if set(x).isdisjoint(['-best_hit_overhang', '-best_hit_score_edge']) else False),
    'db_hard_mask':        (int, None, lambda x: True if set(x).isdisjoint(['-db_soft_mask', '-subject', '-subject_loc']) else False),
    'db_soft_mask':        (int, None, lambda x: True if set(x).isdisjoint(['-db_hard_mask', '-subject', '-subject_loc']) else False),
    'dbsize':        (int, None, None),
    'entrez_query':        (str, None, lambda x: False if set(x).isdisjoint(['-remote']) else True),
    'evalue':        (float, None),
    'export_search_strategy':        (str, None, lambda x: True if set(x).isdisjoint(['-import_search_strategy']) else False),
    'gilist':        (str, lambda x: os.path.isfile(x), lambda x: True if set(x).isdisjoint(['-negative_gilist', '-seqidlist', '-remote', '-subject']) else False),
    'import_search_strategy':        (str, None, lambda x: True if set(x).isdisjoint(['-export_search_strategy']) else False),
    'lcase_masking':        (bool, None, None), #flag
    'line_length':        (int, lambda x: x>=1, None),
    'max_hsps':        (int, lambda x: x>=1, None),
    'max_target_seqs':        (int, lambda x: x>=1, lambda x: True if set(x).isdisjoint(['-num_descriptions', '-num_alignments']) else False), 
    'negative_gilist':        (str, lambda x: os.path.isfile(x), lambda x: True if set(x).isdisjoint(['-gilist', '-seqidlist', '-remote', '-subject', '-subject_loc']) else False),
    'num_alignments':        (int, lambda x: x>=0, lambda x: True if set(x).isdisjoint(['-max_target_seqs']) else False),
    'num_descriptions':        (int, lambda x: x>=0, lambda x: True if set(x).isdisjoint(['-max_target_seqs']) else False),
    'num_threads':        (int, lambda x: x>=1, lambda x: True if set(x).isdisjoint(['-remote']) else False),
    'parse_deflines':        (bool, None, None), #flag
    'qcov_hsp_perc':        (float, lambda x: x>=0.0 and x<=100.0, None),
    'remote':        (str, None, lambda x: True if set(x).isdisjoint(['-gilist', '-seqidlist', '-negative_gilist', '-subject_loc']) else False),
    'searchsp':        (int, lambda x: x>=0, None),
    'seqidlist':        (str, None, lambda x: True if set(x).isdisjoint(['-gilist', '-negative_gilist', '-remote', '-subject']) else False),
    'show_gis':        (bool, None, None),
    'soft_masking':        (bool, None, None), #bool
    'subject':        (str, lambda x: os.path.isfile(x),  lambda x: True if set(x).isdisjoint(['-db', '-gilist', '-seqidlist', '-negative_gilist']) else False),
    'subject_loc':        (str, lambda x: True if re.search('[0-9]*\-[0-9]', x) else False, lambda x: True if set(x).isdisjoint(['-db', '-gilist', '-seqidlist', '-negative_gilist']) else False), #not usable with ['subject', 'subject_loc']
    'sum_stats':        (bool, None, None), #bool
    'window_size':        (int, lambda x: x>=0, None),
    'word_size':        (int, {'blastn': (lambda x: x>=4), 'others': (lambda x: x>=2)} , None),
    'xdrop_ungap':         (float, None, None),
    'query_loc':        (str, lambda x: True if re.search('[0-9]*\-[0-9]', x) else False, None),
    #'db':        (str, None,  lambda x: True if set(x).isdisjoint(['-subject', '-subject_loc']) else False),
    #'html':        (bool, None, None),
    #'out':        (str, None, None),
    #'outfmt':        (int, [1,2,3,4,5,6,7,8,9,10,11] , None),
    #'outfmt_5-11':     (str, ['qseqid', 'qgi', 'qacc', 'sseqid', 'sallseqid', 'sgi', 'sallgi', 'sacc', 'sallacc', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq', 'evalue', 'bitscore', 'score', 'length', 'pident', 'nident', 'mismatch', 'positive', 'gapopen', 'gaps', 'ppos', 'frames', 'qframe', 'sframe', 'btop', 'staxids', 'sscinames', 'scomnames', 'sblastnames', 'sskingdoms', 'stitle', 'salltitles', 'sstrand', 'qcovs', 'qcovhsp', 'qcovus'], None, None),
    #'query':        (str, None, None),
    #'h':        (None,None, None),
    #'help':        (None,None, None),
    #'version':        (None,None, None)

}

#Contains all arguments the blastn tools
_blastn_args = {
    'strand':        (str, ['both','plus','minus'], None),
    'gapopen':        (int, None, None),
    'gapextend':        (int, None, None),
    'penalty':        (int, None, None),
    'reward':        (int, None, None),
    'use_index':        (bool, None, None), #bool
    'index_name':        (str, None, None),
    'dust':        (str, None, None),
    'filtering_db':        (str, None, None),
    'window_masker_taxid':        (int, None, None),
    'window_masker_db':        (str, None, None),
    'perc_identity':        (int, lambda x: x>=0.0 and x<=100.0, None),
    'template_type':        (str, ['coding', 'optimal', 'coding_optimal'], lambda x: False if set(x).isdisjoint(['-template_length']) else True),
    'template_length':        (int, [16, 18, 21], lambda x: False if set(x).isdisjoint(['-template_type']) else True),
    'xdrop_gap':        (float, None, None),
    'xdrop_gap_final':        (float, None, None),
    'no_greedy':        (bool, None, None), #flag
    'min_raw_gapped_score':        (int, None, None),
    'ungapped':        (bool, None, None), #flag
    'off_diagonal_range':        (int, lambda x: x>=0, None)
}

#Contains all arguments the blastp tools
_blastp_args = {
    'gapopen':        (int, None, None),
    'gapextend':        (int, None, None),
    'matrix':        (str, None, None),
    'threshold':        (int, None, None),
    'comp_based_stats':        (str, None, None),
    'seg':        (str, ['yes', 'no', 'window locut hitcut'], None),
    'xdrop_gap':        (float, None, None),
    'xdrop_gap_final':        (float, None, None),
    'use_sw_tback':        (bool, None, None) #flag
}

#Contains all arguments the blastx tools
_blastx_args = {
    'strand':        (str, ['both','plus','minus'], None),
    'query_gencode':        (int, lambda x:((x>=1 and x<=6) or (x>=9 and x<=16) or (x>=21 and x<=25)), None),
    'gapopen':        (int, None, None),
    'gapextend':        (int, None, None),
    'max_intron_length':        (int, lambda x: x>=0, None),
    'matrix':        (str, None, None),
    'threshold':        (float, lambda x: x>=0.0, None),
    'comp_based_stats':        (str, None, None),
    'seg':        (str, ['yes', 'no', 'window locut hitcut'], None),
    'xdrop_gap':        (float, None, None),
    'xdrop_gap_final':        (float, None, None),
    'ungapped':        (bool, None, None), #flag
    'use_sw_tback':        (bool, None, None) #flag
}

#Contains all arguments the tblastn tools
_tblastn_args = {
    'gapopen':        (int, None, None),
    'gapextend':        (int, None, None),
    'db_gencode':        (int, lambda x:((x>=1 and x<=6) or (x>=9 and x<=16) or (x>=21 and x<=25)), None),
    'max_intron_length':        (int, lambda x: x>=0, None),
    'matrix':        (str, None, None),
    'threshold':        (int, lambda x: x>=0, None),
    'comp_based_stats':        (str, None, None),
    'seg':        (str, ['yes', 'no', 'window locut hitcut'], None),
    'xdrop_gap':        (float, None, None),
    'xdrop_gap_final':        (float, None, None),
    'ungapped':        (bool, None, None), #flag
    'use_sw_tback':        (bool, None, None), #flag
    'in_pssm':        (bool, lambda x: os.path.isfile(x), lambda x: True if set(x).isdisjoint(['-remote', '-query', '-query_loc']) else False)

}

#Contains all arguments the tblastx tools
_tblastx_args = {
    'strand':        (str, None, None),
    'query_gencode':        (int, lambda x:((x>=1 and x<=6) or (x>=9 and x<=16) or (x>=21 and x<=25)), None),
    'max_intron_length':        (int, lambda x: x>=0, None),
    'matrix':        (str, None, None),
    'threshold':        (float, lambda x: x>=0, None),
    'db_gencode':        (int, lambda x:((x>=1 and x<=6) or (x>=9 and x<=16) or (x>=21 and x<=25)), None),
    'seg':        (str, None, None)
}

#Contains all arguments the makeblast tools
_mkblstdbargs = {
    'input_type':        (str, ['asn1_bin', 'asn1_txt', 'blastdb', 'fasta'], None),
    'title':        (str, None, None),
    'parse_seqids':        (bool, None, None),
    'hash_index':        (bool, None, None),
    'mask_data':        (str, None, None),
    'mask_id':        (str, None, None),
    'mask_desc':        (str, None, None),
    'gi_mask':        (bool, None, None),
    'gi_mask_name':        (str, None, None),
    'max_file_sz':        (str, None, None),
    'logfile':        (str, None, None),
    'taxid':        (int , lambda x: x>=0, None, None),
    'taxid_map':        (str, None, None)
    #'dbtype':        (str, ['nucl', 'prot'], None),
    #'in': :        (str, None),
    #'out':        (str, None),
}


mkblstdb_intype = {
    'asn1_txt':     ['.asn', '.asn1'],
    #'asn1_bin':     [],
    'fasta': ['.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.frn', '.nt'],
    'blastdb': ['.aso'],
}

class _Blast(_utils._App):
    @_utils._protect_run
    def _check_path(self, path):
        # test the "help" option to ensure that a blast exist
        cmd = (path, '-version')
        try:
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = p.communicate()
            if len(stderr): return stderr
        except OSError, e:
            return e.strerror

        if 'blastn' in stdout:
            mo = re.search('blastn: (\d+)\.(\d+)\.(\d+)', stdout)
        elif 'blastp' in stdout:
            mo = re.search('blastp: (\d+)\.(\d+)\.(\d+)', stdout)
        elif 'blastx' in stdout:
            mo = re.search('blastx: (\d+)\.(\d+)\.(\d+)', stdout)
        elif 'tblastn' in stdout:
            mo = re.search('tblastn: (\d+)\.(\d+)\.(\d+)', stdout)
        elif 'tblastx' in stdout:
            mo = re.search('tblastx: (\d+)\.(\d+)\.(\d+)', stdout)
        else:
            mo = re.search('makeblastdb: (\d+)\.(\d+)\.(\d+)', stdout)

        M = int(mo.group(1))
        m = int(mo.group(2))
        if not mo: return 'cannot read version number from Blast'
        if M >= 2 or m >= 4: return None
        else: return 'Blast version ??? or higher is required'


_app_n = _Blast(key='blastn', default='blastn')
_utils.paths._add(_app_n)
_app_p = _Blast(key='blastp', default='blastp')
_utils.paths._add(_app_p)
_app_x = _Blast(key='blastx', default='blastx')
_utils.paths._add(_app_x)
_app_tx = _Blast(key='tblastx', default='tblastx')
_utils.paths._add(_app_tx)
_app_tn = _Blast(key='tblastn', default='tblastn')
_utils.paths._add(_app_tn)
_app_mkdb = _Blast(key='makeblastdb', default='makeblastdb')
_utils.paths._add(_app_mkdb)


def blast_config(bin_path):
    """
    The method 'blast_config' allow to set the application paths of the differents tools of 'Blast'. With the path of the executable directory of Blast tool.
    :bin_path: path of the executable directory of Blast tool, containing executable files of tools: blastn, blastp, blastx, tblastn, tblastx and makeblastdb
    """
    apps = [ _app_n.key, _app_p.key, _app_x.key, _app_tx.key, _app_tn.key, _app_mkdb.key ]
    if not os.path.exists(bin_path): raise  ValueError, 'The directory {0} doesn`t exist'.format(bin_path)
    for app in apps:
        path = bin_path + app
        _utils.paths[app] = bin_path + app

def makeblastdb(data, out_grp=True, out=None, out_path=None , dbtype='nucl'): #, verbose=False
    """
    The makeblastdb application produces BLAST databases from FASTA,ASN or Blastdb files.
    The database created is stored in a folder according the argument 'out_path'
    :data: file containing sequences
    :out: Database name
    :out_path: Database directory. If 'None', the folder Database is create at the current directory
    :dbtype: type of the sequence in the inputfile. 'Nucl' for nucleic acid and 'prot' for protein
    :return: the databath path

    """
    #Folder creation
    if out is None: out ='tmp_blast_database' 
    if out_path:
        if os.path.exists(out_path):
            folder_name = os.path.join(out_path, out)
        else:
            raise  ValueError, 'The directory {0} doesn`t exist'.format(out_path)
    else:
            folder_name = out

    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    #Check if app path is OK
    path = _app_mkdb.get_path()
    if path is None: raise RuntimeError, 'Makeblastdb program not available -- please configure path'

    #Datebase 
    if not isinstance(data, (_interface.Container, _interface.Align, str)): raise TypeError, 'invalid type for `source` argument'
    if isinstance(data, (_interface.Container, _interface.Align)):
        fname = 'tmp_fasta.fsa'
	mapping = {}
        _utils._write(data, out_grp, fname , mapping)
        if not os.path.isfile(fname): raise  ValueError, 'The file {0} doesn`t exist'.format(fname)
    else:
        if isinstance(data, str) and not (os.path.isfile(data)) : raise ValueError, 'It is not a file' #A modifier
	#Check file type
        filename, file_ext = os.path.splitext(data)
        if any(file_ext in exts for exts in mkblstdb_intype.values()):
	    fname = data
        else:
           raise ValueError, 'The file {0} has a wrong format'.format(data)

    command_line = [path, '-in', fname]

    #Check dbtype
    if dbtype not in ['nucl', 'prot']: raise ValueError, 'invalid option: {0}'.format(dbtype)
    command_line.extend(['-dbtype', dbtype ])

    #Check output
    path_db = os.path.join(folder_name, out)
    command_line.extend(['-out', path_db])
    p = subprocess.Popen(command_line, stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)#"(None if verbose else subprocess.PIPE))
    stdout, stderr = p.communicate('Y\n')

    ext='.nhr' if dbtype == 'nucl' else '.phr'
    if (not os.path.isfile(path_db+ext)) or (not os.path.isfile(path_db+ext)): 
        if stderr is not None:
            lines = stderr.splitlines()
            if len(lines) and 'error' in lines[-1]:
                raise RuntimeError, 'error while running makeblast: {0}'.format(lines[-1].strip())
        raise RuntimeError, 'unknown error while running makeblast (try running in verbose mode)'
    else:
        return os.path.realpath(path_db)



@_utils._protect_run
def blastn(query, database, outgrp=True, task='blastn', **kwargs): #output, verbose=False,
    """
    The blastn application searches a nucleotide query against nucleotide subject sequences or a nucleotide database. 
    :query:  Query file name or a :class:`~._interface.SampleView` instance containing a nucleotide sequence or multiple nucleotide sequences.
    :database: a BLAST nucleotide database name, a :class:`~._interface.Align`, :class:`~._interface.Container` instance or a fasta file with 
               nucleotide sequence.
    :outgrp: boolean indicates if the database created from a :class:`~._interface.Align` or :class:`~._interface.Container`
             instance must include outgroupe.
    :kwargs: other keyword arguments are passed to Blastn. The available options are listed below:
             .. include:: blastn_arguments.txt
    :task:     Four different programs are supported: 
                    `megablast`, for very similar sequences (e.g, sequencing errors),
                    `dc-megablast`, typically used for inter-species comparisons, 
                    `blastn`, the traditional program used for inter-species comparisons, 
                    `blastn-short`, optimized for sequences less than 30 nucleotides.
    :return: a :class:`~._wrappers.Blast` instance
    """
    #Blastn Configuration
    path = _app_n.get_path()
    if path is None: raise RuntimeError, 'Blastn program not available -- please configure path'

    #Tools choice
    tasks_allowed=['blastn', 'blastn-short', 'dc-megablast', 'megablast']
    if not isinstance(task, str): raise TypeError, 'invalid value for argument task: {0}'.format(task)
    if task not in tasks_allowed: raise ValueError, 'invalid value for argument {0}'.format(task)
    else: command_line = [path, '-task', task]

    #Query
    if not isinstance(query, (_interface.SampleView, _interface.Container, _interface.Align, str)): raise TypeError, 'invalid type for `query` argument'
    if isinstance(query, (_interface.SampleView)): 
        tmp_query = 'tmp_query.fsa'
        try: _utils._write_sample(query, tmp_query)
        except (RuntimeError): raise RuntimeError, 'An error occur with the argument query: {0}'.format(database)
        command_line.extend(['-query', tmp_query])
    elif isinstance(query , str) and (os.path.isfile(query)):
	try: _fasta.from_fasta(query, cls=_interface.Container)
        except (ValueError): raise ValueError, 'the file:``{0}`` as query argument is not a fasta file '.format(query)
        else: command_line.extend(['-query', query])
    else: 
        raise ValueError, 'invalid value for argument {0}'.format(query)

    #Datebase
    if (isinstance(database, (_interface.Container, _interface.Align))) or (isinstance(database, str) and (os.path.isfile(database))):
        tmp_db = 'tmp_blast_db'
        try: tmp_db = makeblastdb(database, out_grp=outgrp, out=tmp_db, dbtype='nucl') #, verbose=False
        except (RuntimeError): raise RuntimeError, 'Cannot create a Blast database from the data: {0}'.format(database)
        command_line.extend(['-db', tmp_db])
    elif os.path.isfile(database+'.nhr') and os.path.isfile(database+'.nin') and os.path.isfile(database+'.nsq'):
        command_line.extend(['-db', database])
    else: 
        raise ValueError, 'invalid value for argument {0}'.format(database)

    #TMP out file "XML"  
    command_line.extend(['-out', 'tmp_result.xml'])
    command_line.extend(['-outfmt' , '5'])

    #Check if args is correct and valuetask
    for opt, value in kwargs.iteritems():
        if opt not in _common_args :
            if opt not in _blastn_args: raise ValueError, 'invalid option: {0}'.format(opt)
        
        if opt in _common_args:
            type_ = _common_args[opt][0]
            test = _common_args[opt][1]
        elif opt in _blastn_args:
            type_ = _blastn_args[opt][0]
            test = _blastn_args[opt][1]
        try: value = type_(value)
        except (TypeError, ValueError): raise ValueError, 'invalid value for argument {0}'.format(opt)
        if type_ == str and type(test) is list or type(test) is tuple:
            if value not in test: raise ValueError, 'invalid value for argument {0}'.format(opt)
        elif test is not None and not test(value): raise ValueError, 'invalid value for argument {0}'.format(opt)

        if opt == 'show_gis':
            if value == True: command_line.append('-show_gis')
        elif opt == 'parse_deflines':
            if value == True: command_line.append('-parse_deflines')
        elif opt == 'remote':
            if value == True: command_line.append('-remote')
        elif opt == 'use_index':
            if value == True:
                command_line.append('-use_index=true')
            else:
                command_line.append('-use_index=false')
        elif opt == 'soft_masking':
            if value == True:
                command_line.append('-soft_masking=true')
            else:
                command_line.append('-soft_masking=false')
        elif opt == 'lcase_masking':
            if value == True: command_line.append('-lcase_masking')
        elif opt == 'sum_stats':
            if value == True:
                command_line.append('-sum_stats=true')
            else:
                command_line.append('-sum_stats=false')
        elif opt == 'no_greedy':
            if value == True: command_line.append('-no_greedy')
        elif opt == 'ungapped':
            if value == True: command_line.append('-ungaped')
        else:
            command_line.extend(['-' + opt, str(value)])
        

    # run the program
    cmd_options = [x.replace('-','') for x in command_line if re.search('^-',x)]
    for opt in cmd_options:
        if opt in _common_args: comp = _common_args[opt][2]
        if opt in _blastn_args: comp = _blastn_args[opt][2]
        if opt in ['query','db','task','outfmt','out'] : comp = None
        if comp is not None and not comp(command_line): raise ValueError, 'Incompatibility of the option {0}'.format(opt)

    p = subprocess.Popen(command_line, stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)  ##(None if verbose else subprocess.PIPE))
    stdout, stderr = p.communicate('Y\n')

    # check error
    if not os.path.isfile('tmp_result.xml') or os.path.getsize('tmp_result.xml') == 0:
        if stderr is not None:
            lines = stderr.splitlines()
            if len(lines) and 'Error' in lines[-1]:
                raise RuntimeError, 'error while running Blast: {0}'.format(lines[-1].strip())
        raise RuntimeError, 'unknown error while running Blast (try running in verbose mode)'
    else:   
        Blast = get_data_xml('tmp_result.xml')
        return Blast


@_utils._protect_run
def blastp(query, database, task='blastp', outgrp=False, **kwargs): 
    """
    The blastp application searches a protein sequence against protein subject sequences or a protein database. 
    :query:  Query file name or a :class:`~._interface.SampleView` instance containing a protein sequence or multiple protein sequences.
    :database: a BLAST protein database name, a :class:`~._interface.Align`, :class:`~._interface.Container` instance or a fasta file with 
               protein sequence.
    :outgrp:  boolean indicates if the database created from a :class:`~._interface.Align` or :class:`~._interface.Container`
              instance must include outgroupe.
    :kwargs:  other keyword arguments are passed to Blastn. The available options are listed below:
              .. include:: blastp_arguments.txt
    :task:    Two different programs are supported: 
                   `blastp` (Default), for standard protein-protein comparisons,
                   `blastp-short`, optimized for query sequences shorter than 30 residues.
    :return: a :class:`~._wrappers.Blast` instance
    """ 
    #Check programm path
    path = _app_p.get_path()
    if path is None: raise RuntimeError, 'Blastn program not available -- please configure path'
  
    #Tools choice
    tasks_allowed=['blastp', 'blastp-fast', 'blastp-short']
    if not isinstance(task, str): raise TypeError, 'invalid value for argument task: {0}'.format(task)
    if task not in tasks_allowed: raise ValueError, 'invalid value for argument {0}'.format(task)
    else: command_line = [path, '-task', task]

    #Query
    if not isinstance(query, (_interface.SampleView, str)): raise TypeError, 'invalid type for `query` argument'
    if isinstance(query, (_interface.SampleView)): 
        tmp_query = 'tmp_query.fsa'
        try: _utils._write_sample(query, tmp_query)
        except (RuntimeError): raise RuntimeError, 'An error occur with the argument query: {0}'.format(database)
        command_line.extend(['-query', tmp_query])
    elif isinstance(query , str) and os.path.isfile(query):
        try: _fasta.from_fasta(query, cls=_interface.Container)
        except ValueError: raise ValueError, 'the file:``{0}`` as query argument is not a fasta file '.format(query)
        else: command_line.extend(['-query', query])
    else: 
        raise ValueError, 'invalid value for argument {0}'.format(query)

    #Datebase
    if (isinstance(database, (_interface.Container, _interface.Align))) or (isinstance(database, str) and (os.path.isfile(database))):
        tmp_db = 'tmp_blast_db'
        try: tmp_db = makeblastdb(database, out_grp=outgrp, out=tmp_db, dbtype='prot') #, verbose=False
        except (RuntimeError): raise RuntimeError, 'Cannot create a Blast database from the data: {0}'.format(database)
        command_line.extend(['-db', tmp_db])
    elif os.path.isfile(database+'.phr') and os.path.isfile(database+'.pin') and os.path.isfile(database+'.psq'):
        command_line.extend(['-db', database])
    else: 
        raise ValueError, 'invalid value for argument {0}'.format(database)
    
    #TMP out file "XML"  
    command_line.extend(['-out', 'tmp_result.xml'])
    command_line.extend(['-outfmt' , '5'])

    #Check if args is correct and value
    for opt, value in kwargs.iteritems():
        if opt not in _common_args :
            if opt not in  _blastp_args: raise ValueError, 'invalid option: {0}'.format(opt)
        
        if opt in _common_args:
            type_ = _common_args[opt][0]
            test = _common_args[opt][1]
        elif opt in _blastp_args:
            type_ = _blastp_args[opt][0]
            test = _blastp_args[opt][1]
        try: value = type_(value)
        except (TypeError, ValueError): raise ValueError, 'invalid value for argument {0}'.format(opt)
        if type_ == str and type(test) is list or type(test) is tuple:
            if value not in test: raise ValueError, 'invalid value for argument {0}'.format(opt)
        elif test is not None and not test(value): raise ValueError, 'invalid value for argument {0}'.format(opt)

        if opt == 'show_gis':
            if value == True: command_line.append('-show_gis')
        elif opt == 'soft_masking':
            if value == True:
                command_line.append('-soft_masking=true')
            else:
                command_line.append('-soft_masking=false')
        elif opt == 'lcase_masking':
            if value == True: command_line.append('-lcase_masking')
        elif opt == 'sum_stats':
            if value == True:
                command_line.append('-sum_stats=true')
            else:
                command_line.append('-sum_stats=false')
        elif opt == 'ungapped':
            if value == True: command_line.append('-ungapped')
        elif opt == 'parse_deflines':
            if value == True: command_line.append('-parse_deflines')
        elif opt == 'remote':
            if value == True: command_line.append('-remote')
        elif opt == 'use_sw_tback': 
            if value == True: command_line.append('-use_sw_tback')
        else:
            command_line.extend(['-' + opt, str(value)])

    # run the program
    cmd_options = [x.replace('-','') for x in command_line if re.search('^-',x)]
    for opt in cmd_options:
        if opt in _common_args: comp = _common_args[opt][2]
        if opt in _blastp_args: comp = _blastp_args[opt][2]
        if opt in ['query','db','task','outfmt','out'] : comp = None
        if comp is not None and not comp(command_line): raise ValueError, 'Incompatibility of the option {0}'.format(opt)

    p = subprocess.Popen(command_line, stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)  ##(None if verbose else subprocess.PIPE))
    stdout, stderr = p.communicate('Y\n')


    # check error
    if not os.path.isfile('tmp_result.xml') or os.path.getsize('tmp_result.xml') == 0:
        if stderr is not None:
            lines = stderr.splitlines()
            if len(lines) and 'Error' in lines[-1]:
                raise RuntimeError, 'error while running Blast: {0}'.format(lines[-1].strip())
        raise RuntimeError, 'unknown error while running Blast (try running in verbose mode)'
    else:   
        Blast = get_data_xml('tmp_result.xml')
        return Blast

@_utils._protect_run
def blastx(query, database ,outgrp=True, task='blastx', **kwargs):
    """
    The blastx application translates a nucleotide query and searches it against protein subject sequences or a protein database
    :query:  Query file name or a :class:`~._interface.SampleView` instance containing a nucleotide sequence or multiple nucleotide sequence.
    :database: a BLAST protein database name, a :class:`~._interface.Align`, :class:`~._interface.Container` instance or a fasta file with 
               protein sequence.
    :database: BLAST database name, a :class:`~._interface.Align` or :class:`~._interface.Container` instance.
    :outgrp: boolean indicates if the database created from a :class:`~._interface.Align` or :class:`~._interface.Container`
             instance must include outgroupe.
    :kwargs: other keyword arguments are passed to Blastn. The available options are listed below:
             .. include:: tblastx_arguments.txt
    :task:    Two different programs are supported: 
                   `blastx` (Default) for standard comparisons
                   `blastx-fast` optimized for query sequences shorter than 30 residues
    :return: a :class:`~._wrappers.Blast` instance
    """

    #Check programm path
    path = _app_x.get_path()
    if path is None: raise RuntimeError, 'Blastx program not available -- please configure path'

    #Tools choice
    tasks_allowed=['blastx', 'blastx-fast']
    if not isinstance(task, str): raise TypeError, 'invalid value for argument task: {0}'.format(task)
    if task not in tasks_allowed: raise ValueError, 'invalid value for argument {0}'.format(task)
    else: command_line = [path, '-task', task]

    #Query
    if not isinstance(query, (_interface.SampleView, str)): raise TypeError, 'invalid type for `query` argument'
    if isinstance(query, (_interface.SampleView)): 
        tmp_query = 'tmp_query.fsa'
        try: _utils._write_sample(query, tmp_query)
        except (RuntimeError): raise RuntimeError, 'An error occur with the argument query: {0}'.format(database)
        command_line.extend(['-query', tmp_query])
    elif isinstance(query , str) and (os.path.isfile(query)):
	try: _fasta.from_fasta(query, cls=_interface.Container)
        except (ValueError): raise ValueError, 'the file:``{0}`` as query argument is not a fasta file '.format(query)
        else: command_line.extend(['-query', query])
    else: 
        raise ValueError, 'invalid value for argument {0}'.format(query)

    #Datebase
    if (isinstance(database, (_interface.Container, _interface.Align))) or (isinstance(database, str) and (os.path.isfile(database))):
        tmp_db = 'tmp_blast_db'
        try: tmp_db = makeblastdb(database, out_grp=outgrp, out=tmp_db, dbtype='prot') #, verbose=False
        except (RuntimeError): raise RuntimeError, 'Cannot create a Blast database from the data: {0}'.format(database)
        command_line.extend(['-db', tmp_db])
    elif os.path.isfile(database+'.phr') and os.path.isfile(database+'.pin') and os.path.isfile(database+'.psq'):
        command_line.extend(['-db', database])
    else: 
        raise ValueError, 'invalid value for argument {0}'.format(database)
      
    #TMP out file "XML"  
    command_line.extend(['-out', 'tmp_result.xml'])
    command_line.extend(['-outfmt' , '5'])

    #Check if args is correct and value
    for opt, value in kwargs.iteritems():
        if opt not in _common_args :
            if opt not in _blastx_args: raise ValueError, 'invalid option: {0}'.format(opt)
        
        if opt in _common_args: type_, test = _common_args[opt]
        elif opt in _blastx_args: type_, test = _blastx_args[opt]
        try: value = type_(value)
        except (TypeError, ValueError): raise ValueError, 'invalid value for argument {0}'.format(opt)
        if type_ == str and type(test) is list or type(test) is tuple:
            if value not in test: raise ValueError, 'invalid value for argument {0}'.format(opt)
        elif test is not None and not test(value): raise ValueError, 'invalid value for argument {0}'.format(opt)

        if opt == 'show_gis':
            if value == True: command_line.append('-show_gis')
        elif opt == 'soft_masking':
            if value == True:
                command_line.append('-soft_masking=true')
            else:
                command_line.append('-soft_masking=false')
        elif opt == 'lcase_masking':
            if value == True: command_line.append('-lcase_masking')
        elif opt == 'sum_stats':
            if value == True:
                command_line.append('-sum_stats=true')
            else:
                command_line.append('-sum_stats=false')
        elif opt == 'ungapped':
            if value == True: command_line.append('-ungapped')
        elif opt == 'parse_deflines':
            if value == True: command_line.append('-parse_deflines')
        elif opt == 'remote':
            if value == True: command_line.append('-remote')
        elif opt == 'use_sw_tback': 
            if value == True: command_line.append('-use_sw_tback')
        else:
            command_line.extend(['-' + opt, str(value)])

    # run the program
    cmd_options = [x.replace('-','') for x in command_line if re.search('^-',x)]
    for opt in cmd_options:
        if opt in _common_args: comp = _common_args[opt][2]
        if opt in _blastn_args: comp = _blastx_args[opt][2]
        if opt in ['query','db','task','outfmt','out'] : comp = None
        if comp is not None and not comp(command_line): raise ValueError, 'Incompatibility of the option {0}'.format(opt)

    p = subprocess.Popen(command_line, stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)  ##(None if verbose else subprocess.PIPE))
    stdout, stderr = p.communicate('Y\n')


    # check error
    if not os.path.isfile('tmp_result.xml') or os.path.getsize('tmp_result.xml') == 0:
        if stderr is not None:
            lines = stderr.splitlines()
            if len(lines) and 'Error' in lines[-1]:
                raise RuntimeError, 'error while running Blast: {0}'.format(lines[-1].strip())
        raise RuntimeError, 'unknown error while running Blast (try running in verbose mode)'
    else:   
        Blast = get_data_xml('tmp_result.xml')
        return Blast


@_utils._protect_run
def tblastn(query, database, outgrp=True, task='tblastn', **kwargs):
    """
    The tblastn application searches a protein query against nucleotide subject sequences or a nucleotide database translated at search time.
    :query:    a fasta file name or a :class:`~._interface.SampleView` instance containing a protein sequence or multiple protein sequence as query.
    :database: a BLAST nucleotide database name, a :class:`~._interface.Align`, :class:`~._interface.Container` instance or a fasta file with 
               nucleotide sequence.
    :outgrp:   boolean indicates if the database created from a :class:`~._interface.Align` or :class:`~._interface.Container`
               instance must include outgroupe.
    :task:    Two different programs are supported:
                   `tblastn`(Default)
                   `tblastn-fast`
    :kwargs:   other keyword arguments are passed to Blastn. The available options are listed below:
               .. include:: tblastn_arguments.txt
    :return: a :class:`~._wrappers.Blast` instance
    """
    path = _app_tn.get_path()
    if path is None: raise RuntimeError, 'tBlastn program not available -- please configure path'

    #Tools choice
    tasks_allowed=['tblastn', 'tblastn-fast']
    if not isinstance(task, str): raise TypeError, 'invalid value for argument task: {0}'.format(task)
    if task not in tasks_allowed: raise ValueError, 'invalid value for argument {0}'.format(task)
    else: command_line = [path, '-task', task]

    #Query
    if not isinstance(query, (_interface.SampleView, str)): raise TypeError, 'invalid type for `query` argument'
    if isinstance(query, (_interface.SampleView)): 
        tmp_query = 'tmp_query.fsa'
        try: _utils._write_sample(query, tmp_query)
        except (RuntimeError): raise RuntimeError, 'An error occur with the argument query: {0}'.format(database)
        command_line.extend(['-query', tmp_query])
    elif isinstance(query , str) and (os.path.isfile(query)):
	try: _fasta.from_fasta(query, cls=_interface.Container)
        except (ValueError): raise ValueError, 'the file:``{0}`` as query argument is not a fasta file '.format(query)
        else: command_line.extend(['-query', query])            
    else: 
        raise ValueError, 'invalid value for argument {0}'.format(query)

    #Datebase
    if (isinstance(database, (_interface.Container, _interface.Align))) or (isinstance(database, str) and (os.path.isfile(database))):
        tmp_db = 'tmp_blast_db'
        try: tmp_db = makeblastdb(database, out_grp=outgrp, out=tmp_db, dbtype='nucl') #, verbose=False
        except (RuntimeError): raise RuntimeError, 'Cannot create a Blast database from the data: {0}'.format(database)
        command_line.extend(['-db', tmp_db])
    elif os.path.isfile(database+'.nhr') and os.path.isfile(database+'.nin') and os.path.isfile(database+'.nsq'):
        command_line.extend(['-db', database])
    else: 
        raise ValueError, 'invalid value for argument {0}'.format(database)
    
    #TMP out file "XML"  
    command_line.extend(['-out', 'tmp_result.xml'])
    command_line.extend(['-outfmt' , '5'])

    #Check if args is correct and value
    for opt, value in kwargs.iteritems():
        if opt not in _common_args :
            if opt not in _tblastn_args: raise ValueError, 'invalid option: {0}'.format(opt)
        
        if opt in _common_args:
            type_ = _common_args[opt][0]
            test = _common_args[opt][1]
        elif opt in _tblastn_args:
            type_ = _tblastn_args[opt][0]
            test = _tblastn_args[opt][1]
        try: value = type_(value)
        except (TypeError, ValueError): raise ValueError, 'invalid value for argument {0}'.format(opt)
        if type_ == str and type(test) is list or type(test) is tuple:
            if value not in test: raise ValueError, 'invalid value for argument {0}'.format(opt)
        elif test is not None and not test(value): raise ValueError, 'invalid value for argument {0}'.format(opt)

        if opt == 'show_gis':
            if value == True: command_line.append('-show_gis')
        elif opt == 'soft_masking':
            if value == True:
                command_line.append('-soft_masking=true')
            else:
                command_line.append('-soft_masking=false')
        elif opt == 'lcase_masking':
            if value == True: command_line.append('-lcase_masking')
        elif opt == 'sum_stats':
            if value == True:
                command_line.append('-sum_stats=true')
            else:
                command_line.append('-sum_stats=false')
        elif opt == 'ungapped':
            if value == True: command_line.append('-ungapped')
        elif opt == 'parse_deflines':
            if value == True: command_line.append('-parse_deflines')
        elif opt == 'remote':
            if value == True: command_line.append('-remote')
        elif opt == 'use_sw_tback': 
            if value == True: command_line.append('-use_sw_tback')
        else:
            command_line.extend(['-' + opt, str(value)])

    # run the program
    cmd_options = [x.replace('-','') for x in command_line if re.search('^-',x)]
    for opt in cmd_options:
        if opt in _common_args: comp = _common_args[opt][2]
        if opt in _blastn_args: comp = _tblastn_args[opt][2]
        if opt in ['query','db','task','outfmt','out'] : comp = None
        if comp is not None and not comp(command_line): raise ValueError, 'Incompatibility of the option {0}'.format(opt)

    p = subprocess.Popen(command_line, stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)  ##(None if verbose else subprocess.PIPE))
    stdout, stderr = p.communicate('Y\n')

    # check error
    if not os.path.isfile('tmp_result.xml') or os.path.getsize('tmp_result.xml') == 0:
        if stderr is not None:
            lines = stderr.splitlines()
            if len(lines) and 'Error' in lines[-1]:
                raise RuntimeError, 'error while running Blast: {0}'.format(lines[-1].strip())
        raise RuntimeError, 'unknown error while running Blast (try running in verbose mode)'
    else:   
        Blast = get_data_xml('tmp_result.xml')
        return Blast


@_utils._protect_run
def tblastx(query, database, outgrp=True, **kwargs):
    """
    The tblastx application searches a translated nucleotide query against translated nucleotide subject sequences or a 
    translated nucleotide database.
    :query:  Query file name or a :class:`~._interface.SampleView` instance containing a nucleotide(translated) sequence or multiple 
             nucleotide(translated) sequence as query.
    :database: a BLAST nucleotide(translated) database name, a :class:`~._interface.Align`, :class:`~._interface.Container` instance 
               or a fasta file with nucleotide(translated) sequence.
    :outgrp: boolean indicates if the database created from a :class:`~._interface.Align` or :class:`~._interface.Container`
             instance must include outgroupe.
    :kwargs: other keyword arguments are passed to Blastn. The available options are listed below:
             .. include:: tblastx_arguments.txt
    :return: a :class:`~._wrappers.Blast` instance
    """
    path = _app_tx.get_path()
    command_line=[path]
    if path is None: raise RuntimeError, 'tBlastx program not available -- please configure path'

    #Query
    if not isinstance(query, (_interface.SampleView, str)): raise TypeError, 'invalid type for `query` argument'
    if isinstance(query, (_interface.SampleView)): 
        tmp_query = 'tmp_query.fsa'
        try: _utils._write_sample(query, tmp_query)
        except (RuntimeError): raise RuntimeError, 'An error occur with the argument query: {0}'.format(database)
        command_line.extend(['-query', tmp_query])
    elif isinstance(query , str) and (os.path.isfile(query)):
	try: _fasta.from_fasta(query, cls=_interface.Container)
        except (ValueError): raise ValueError, 'the file:``{0}`` as query argument is not a fasta file '.format(query)
        else: command_line.extend(['-query', query])
    else: 
        raise ValueError, 'invalid value for argument {0}'.format(query)

    #Datebase
    if (isinstance(database, (_interface.Container, _interface.Align))) or (isinstance(database, str) and (os.path.isfile(database))):
        tmp_db = 'tmp_blast_db'
        try: tmp_db = makeblastdb(database, out_grp=outgrp, out=tmp_db, dbtype='nucl') #, verbose=False
        except (RuntimeError): raise RuntimeError, 'Cannot create a Blast database from the data: {0}'.format(database)
        command_line.extend(['-db', tmp_db])
    elif os.path.isfile(database+'.nhr') and os.path.isfile(database+'.nin') and os.path.isfile(database+'.nsq'):
        command_line.extend(['-db', database])
    else: 
        raise ValueError, 'invalid value for argument {0}'.format(task)
        
    #TMP out file "XML" 
    command_line.extend(['-out', 'tmp_result.xml'])
    command_line.extend(['-outfmt' , '5'])

    #Check if args is correct and value
    for opt, value in kwargs.iteritems():
        if opt not in _common_args :
            if opt not in _tblastx_args: raise ValueError, 'invalid option: {0}'.format(opt)
        
        if opt in _common_args:
            type_ = _common_args[opt][0]
            test = _common_args[opt][1]
        elif opt in _tblastx_args:
            type_ = _tblastx_args[opt][0]
            test = _tblastx_args[opt][1]
        try: value = type_(value)
        except (TypeError, ValueError): raise ValueError, 'invalid value for argument {0}'.format(opt)
        if type_ == str and type(test) is list or type(test) is tuple:
            if value not in test: raise ValueError, 'invalid value for argument {0}'.format(opt)
        elif test is not None and not test(value): raise ValueError, 'invalid value for argument {0}'.format(opt)

        if opt == 'show_gis':
            if value == True: command_line.append('-show_gis')
        elif opt == 'soft_masking':
            if value == True:
                command_line.append('-soft_masking=true')
            else:
                command_line.append('-soft_masking=false')
        elif opt == 'lcase_masking':
            if value == True: command_line.append('-lcase_masking')
        elif opt == 'sum_stats':
            if value == True:
                command_line.append('-sum_stats=true')
            else:
                command_line.append('-sum_stats=false')
        elif opt == 'parse_deflines':
            if value == True: command_line.append('-parse_deflines')
        elif opt == 'remote':
            if value == True: command_line.append('-remote')
        else:
            command_line.extend(['-' + opt, str(value)])

    # run the program
    cmd_options = [x.replace('-','') for x in command_line if re.search('^-',x)]
    for opt in cmd_options:
        if opt in _common_args: comp = _common_args[opt][2]
        if opt in _blastn_args: comp = _tblastx_args[opt][2]
        if opt in ['query','db','outfmt','out'] : comp = None
        if comp is not None and not comp(command_line): raise ValueError, 'Incompatibility of the option {0}'.format(opt)

    p = subprocess.Popen(command_line, stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)  ##(None if verbose else subprocess.PIPE))
    stdout, stderr = p.communicate('Y\n')

    # check error
    if not os.path.isfile('tmp_result.xml') or os.path.getsize('tmp_result.xml') == 0:
        if stderr is not None:
            lines = stderr.splitlines()
            if len(lines) and 'Error' in lines[-1]:
                raise RuntimeError, 'error while running Blast: {0}'.format(lines[-1].strip())
        raise RuntimeError, 'unknown error while running Blast (try running in verbose mode)'
    else:   
        Blast = get_data_xml('tmp_result.xml')
        return Blast








