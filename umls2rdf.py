#! /usr/bin/env python

import logging
import sys
import os
from os import path
import urllib
from string import Template
import collections
import optparse
from codecs import open

ENC = 'utf-8'

FORMAT = '%(asctime)-15s %(message)s'
logging.basicConfig(format=FORMAT)
LOG = logging.getLogger(name='umls2rdf')
LOG.setLevel(logging.DEBUG)

KNOWN_ONT_CODES = set(['AIR', 'CST', 'CSP', 'HCPCS', 'ICD10PCS', 
        'ICD10CM', 'ICPC2P', 'ICD9CM', 'LNC', 'MDR', 'MDDB', 'MSH',
        'MEDLINEPLUS', 'MTHCH', 'NDDF', 'NDFRT', 'OMIM', 'PDQ', 
        'RCD', 'RXNORM', 'SNOMEDCT', 'VANDF', 'WHO', 'HL7', 
        'ICD10', 'CPT', 'ICPC']) 

UMLS_BASE_URL = "http://bioportal.bioontology.org/ontologies/umls/"

PREFIXES = """
@prefix skos: <http://www.w3.org/2004/02/skos/core#> .
@prefix owl:  <http://www.w3.org/2002/07/owl#> .
@prefix rdfs:  <http://www.w3.org/2000/01/rdf-schema#> .
@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .
@prefix umls: <%(umls)s> .
@prefix umlssty: <%(umls)ssty/> .
@base <http://purl.bioontology.org/ontology/> .


""" % {'umls': UMLS_BASE_URL}

ONTOLOGY_HEADER = Template("""
<$uri>
    a owl:Ontology ;
    rdfs:comment "$comment" ;
    rdfs:label "$label" ;
    owl:imports <http://www.w3.org/2004/02/skos/core> ;
    owl:versionInfo "$versioninfo" .

""")

STY_URL = "umlssty:"
HAS_STY = "umls:hasSTY"
HAS_AUI = "umls:aui"
HAS_CUI = "umls:cui"
HAS_TUI = "umls:tui"

# http://www.nlm.nih.gov/research/umls/sourcereleasedocs/current/SNOMEDCT/relationships.html

class IDXS(object):

    class MRCONSO(object):
        CODE = 13
        AUI = 7
        STR = 14
        STT = 4
        SCUI = 9
        ISPREF = 6
        TTY = 12
        TS = 2
        CUI = 0
        LAT = 1
        SAB = 11

    class MRREL(object):
        AUI1 = 1
        AUI2 = 5
        CUI1 = 0
        CUI2 = 4
        REL = 3
        RELA = 7
        SAB = 10

    class MRDEF(object):
        AUI = 1
        DEF = 5
        CUI = 0
        SAB = 4

    class MRSAT(object):
        CUI = 0
        CODE = 5
        ATV = 10
        ATN = 8
        SAB = 9

    class MRDOC(object):
        VALUE = 1
        TYPE = 2
        DESC = 3

    class MRRANK(object):
        TTY = 2
        RANK = 0
        SAB = 1

    class MRSTY(object):
        CUI = 0
        TUI = 1
        STN = 2
        STY = 3


def get_umls_url(code):
    return "%s" % code

def flatten(matrix):
    return reduce(lambda x, y: x + y, matrix)

def escape(string):
    return string.replace("\\", "\\\\").replace('"', '\\"')

def get_url_term(ns, code):
    sep = '' if ns[-1] in ('/', ':') else '/'
    ret = "<%s%s%s>"% (ns, sep, urllib.quote(code))
    return ret.replace("%20", "+") 

def get_rel_fragment(rel):
    return rel[IDXS.MRREL.RELA] if rel[IDXS.MRREL.RELA] else rel[IDXS.MRREL.REL]


# NOTE: See UmlsOntology.terms() for the reason these functions use -1 and -2
# indices to obtain the source and target codes, respectively.
def get_rel_code_source(rel, on_cuis):
    return rel[-1] if not on_cuis else rel[IDXS.MRREL.CUI2]
def get_rel_code_target(rel, on_cuis):
    return rel[-2] if not on_cuis else rel[IDXS.MRREL.CUI1]

def get_code(reg, load_on_cuis):
    if load_on_cuis:
        return reg[IDXS.MRCONSO.CUI]
    if reg[IDXS.MRCONSO.CODE]:
        return reg[IDXS.MRCONSO.CODE]
    raise AttributeError, "No code on reg [%s]"%("|".join(reg))

def generate_semantic_types(data_root, url, fileout):
    hierarchy = collections.defaultdict(lambda : list()) 
    all_nodes = list()
    mrsty = UmlsTable("MRSTY", data_root, columns=['TUI', 'STN', 'STY'])
    ont = list()
    ont.append(PREFIXES)
    types_seen = set()
    for stt in mrsty.scan():
        if tuple(stt) in types_seen:
            continue
        types_seen.add(tuple(stt))
        hierarchy[stt[1]].append(stt[0])
        sty_term = """%s a owl:Class ;
\tskos:notation "%s"^^xsd:string ;
\tskos:prefLabel "%s"@en .
""" % (url+stt[0], stt[0], stt[2])
        ont.append(sty_term)
        all_nodes.append(stt)
    
    for node in all_nodes:
        parent = ".".join(node[1].split(".")[0:-1])
        rdfs_subclasses = ["<%s> rdfs:subClassOf %s ." % 
            (url+node[0], url+x) for x in hierarchy[parent]]
        for sc in rdfs_subclasses:
            ont.append(sc)
    data_ont_ttl = "\n".join(ont)
    with open(fileout, "w", encoding=ENC) as fout:
        fout.write(data_ont_ttl)
        fout.write("\n")
        fout.close()

_CACHED_TABLES = {}

class UmlsTable(object):
    def __init__(self, table_name, data_root, columns=None,
            cacheable_sabs=None):
        self.table_name = table_name
        self.data_root = data_root
        self._indexes = getattr(IDXS, table_name)
        self.columns = columns
        self.cacheable_sabs = cacheable_sabs
        self._cached_rows = None

    def _col_idx(self, colname):
        return getattr(self._indexes, colname)

    def scan(self, filt=None):
        for row in self._get_rows():
            if filt and any(row[self._col_idx(col)] != tgt
                    for col, tgt in filt.iteritems()):
                continue
            if self.columns:
                yield [row[self._col_idx(cn)] for cn in self.columns]
            else:
                yield row

    def _get_rows(self):
        if self._cached_rows is not None:
            for row in self._cached_rows:
                yield row
        else:
            caching = self.cacheable_sabs is not None
            if caching:
                self._cached_rows = []
            rrf_path = path.join(self.data_root, "%s.RRF" % self.table_name)
            with open(rrf_path, encoding=ENC) as f:
                for line in f:
                    elems = line.rstrip('\n').split('|')
                    if caching and elems[self._col_idx('SAB')] in self.cacheable_sabs:
                        self._cached_rows.append(elems)
                    yield elems

    @staticmethod
    def get(table_name, data_root, columns=None, cacheable_sabs=None):
        if cacheable_sabs is None:
            return UmlsTable(table_name, data_root, columns)
        key = (table_name, data_root, columns, cacheable_sabs)
        if key not in _CACHED_TABLES:
            _CACHED_TABLES[key] = UmlsTable(table_name, data_root, columns, cacheable_sabs)
        return _CACHED_TABLES[key]


class UmlsClass(object):
    def __init__(self, ns, atoms=None, rels=None,
            defs=None, atts=None, rank=None,
            rank_by_tty=None, sty=None,
            sty_by_cui=None, load_on_cuis=False,
            is_root=None, ont_link=True):
        self.ns = ns
        self.atoms = atoms
        self.rels = rels
        self.defs = defs
        self.atts = atts
        self.rank = rank
        self.rank_by_tty = rank_by_tty
        self.sty = sty
        self.sty_by_cui = sty_by_cui
        self.load_on_cuis = load_on_cuis
        self.is_root = is_root
        self.ont_link = ont_link

    def code(self):
        codes = set([get_code(x, self.load_on_cuis) for x in self.atoms])
        if len(codes) != 1:
            raise AttributeError, "Only one code per term."
        #if DEBUG:
            #LOG.debug(self.atoms)
            #LOG.debug(codes)
        return codes.pop()

    def getAltLabels(self, prefLabel):
        #is_pref_atoms =  filter(lambda x: x[IDXS.MRCONSO.ISPREF] == 'Y', self.atoms)
        return set([atom[IDXS.MRCONSO.STR]
                for atom in self.atoms 
                if atom[IDXS.MRCONSO.STR] != prefLabel])
        
    def getPrefLabel(self):
        if self.load_on_cuis:
            if len(self.atoms) == 1:
                return self.atoms[0][IDXS.MRCONSO.STR]

            labels = set([x[IDXS.MRCONSO.STR] for x in self.atoms])
            if len(labels) == 1:
                return labels.pop()

            #if there's only one ISPREF=Y then that one.
            is_pref_atoms =  filter(lambda x: x[IDXS.MRCONSO.ISPREF] == 'Y', self.atoms)
            if len(is_pref_atoms) == 1:
                return is_pref_atoms[0][IDXS.MRCONSO.STR]
            elif len(is_pref_atoms) > 1:
                is_pref_atoms =  filter(lambda x: x[IDXS.MRCONSO.STT] == 'PF', is_pref_atoms)
                if len(is_pref_atoms) > 0:
                    return is_pref_atoms[0][IDXS.MRCONSO.STR]
            is_pref_atoms =  filter(lambda x: x[IDXS.MRCONSO.STT] == 'PF', self.atoms)
            if len(is_pref_atoms) == 1:
                return is_pref_atoms[0][IDXS.MRCONSO.STR]
            return self.atoms[0][IDXS.MRCONSO.STR]
        else:
            #if ISPREF=Y is not 1 then we look into MRRANK.
            if len(self.rank) > 0:
                sort_key = \
                lambda x: int(self.rank[self.rank_by_tty[x[IDXS.MRCONSO.TTY]][0]][IDXS.MRRANK.RANK])
                mmrank_sorted_atoms = sorted(self.atoms, key=sort_key, reverse=True)
                return mmrank_sorted_atoms[0][IDXS.MRCONSO.STR]
            #there is no rank to use
            else:
                pref_atom = filter(lambda x: 'P' in x[IDXS.MRCONSO.TTY], self.atoms)
                if len(pref_atom) == 1:
                    return pref_atom[0][IDXS.MRCONSO.STR]
            raise AttributeError, "Unable to select pref label"
    
    def getURLTerm(self, code):
        return get_url_term(self.ns, code)
    
    def toRDF(self, fmt="Turtle", hierarchy=True):
        if not fmt == "Turtle":
            raise AttributeError, "Only fmt='Turtle' is currently supported"
        term_code = self.code()
        url_term = self.getURLTerm(term_code)
        prefLabel = self.getPrefLabel()
        altLabels = self.getAltLabels(prefLabel)
        rdf_term = """%s a owl:Class ;
\tskos:prefLabel \"\"\"%s\"\"\"@en ;
\tskos:notation \"\"\"%s\"\"\"^^xsd:string ;
""" % (url_term, escape(prefLabel), escape(term_code))
        if len(altLabels) > 0:
            rdf_term += """\tskos:altLabel %s ;
""" % (" , ".join('\"\"\"%s\"\"\"@en' % escape(al) for al in altLabels))
        if self.is_root: 
            rdf_term += '\tumls:isRoot "true"^^xsd:boolean ;\n'
            # TODO: Discuss adding this subclass relation.
            #rdf_term += '\trdfs:subClassOf owl:Thing ;\n'

        if len(self.defs) > 0:
            rdf_term += """\tskos:definition %s ;
""" % (" , ".join('\"\"\"%s\"\"\"@en' % escape(d[IDXS.MRDEF.DEF]) for d in self.defs))

        for rel in self.rels:
            source_code = get_rel_code_source(rel, self.load_on_cuis)
            target_code = get_rel_code_target(rel, self.load_on_cuis)
            if source_code != term_code:
                raise AttributeError, "Inconsistent code in rel"
            # Map child relations to rdf:subClassOf (skip parent relations).
            if rel[IDXS.MRREL.REL] == 'PAR':
                continue
            if rel[IDXS.MRREL.REL] == 'CHD' and hierarchy:
                o = self.getURLTerm(target_code)
                rdf_term += "\trdfs:subClassOf %s ;\n" % (o, )
            else:
                p = self.getURLTerm(get_rel_fragment(rel))
                o = self.getURLTerm(target_code)
                rdf_term += "\t%s %s ;\n" % (p, o)

        for att in self.atts:
            atn = att[IDXS.MRSAT.ATN]
            atv = att[IDXS.MRSAT.ATV]
            if atn == 'AQ':
                # Skip all these values (they are replicated in MRREL for
                # SNOMEDCT, unknown relationship for MSH).
                #if DEBUG:
                #  LOG.debug("att: %s\n" % str(att))
                #  sys.stderr.flush()
                continue
            rdf_term += "\t%s \"\"\"%s\"\"\"^^xsd:string ;\n" % (self.getURLTerm(atn), escape(atv))
        if self.ont_link:
            rdf_term += "\trdfs:isDefinedBy <%s> ;\n" % (self.ns,)

        #auis = set([x[IDXS.MRCONSO.AUI] for x in self.atoms])
        cuis = set([x[IDXS.MRCONSO.CUI] for x in self.atoms])
        sty_recs = flatten([indexes for indexes in [self.sty_by_cui[cui] for cui in cuis]])
        types = [self.sty[index][IDXS.MRSTY.TUI] for index in sty_recs]

        #for t in auis:
        #    rdf_term += """\t%s \"\"\"%s\"\"\"^^xsd:string ;\n"""%(HAS_AUI, t)
        for t in cuis:
            rdf_term += """\t%s \"\"\"%s\"\"\"^^xsd:string ;\n""" % (HAS_CUI, t)
        for t in set(types):
            rdf_term += """\t%s \"\"\"%s\"\"\"^^xsd:string ;\n""" % (HAS_TUI, t)
        for t in set(types):
            rdf_term += """\t%s %s ;\n""" % (HAS_STY, STY_URL+t)

        return rdf_term + " .\n\n"



class UmlsAttribute(object):
    def __init__(self, ns, att):
        self.ns = ns
        self.att = att

    def getURLTerm(self, code):
        return get_url_term(self.ns, code)

    def toRDF(self, fmt="Turtle"):
        if not fmt == "Turtle":
            raise AttributeError, "Only fmt='Turtle' is currently supported"
        return """<%s> a owl:DatatypeProperty ;
\trdfs:label \"\"\"%s\"\"\";
\trdfs:comment \"\"\"%s\"\"\" .
\n"""%(self.getURLTerm(self.att[IDXS.MRDOC.VALUE]), escape(self.att[IDXS.MRDOC.VALUE]), escape(self.att[IDXS.MRDOC.DESC]))



class UmlsOntology(object):
    def __init__(self, ont_code, ns, data_root, load_on_cuis=False,
            umls_version="2012AB", store_atts=True, all_ont_codes=None,
            ont_link=True):
        self.loaded = False
        self.ont_code = ont_code
        self.ns = ns
        self.data_root = data_root
        self.load_on_cuis = load_on_cuis
        #self.alt_uri_code = alt_uri_code
        self.atoms = list()
        self.atoms_by_code = collections.defaultdict(list)
        if not self.load_on_cuis:
            self.atoms_by_aui = collections.defaultdict(list)
        self.rels = list()
        self.rels_by_aui_src = collections.defaultdict(list)
        self.defs = list()
        self.defs_by_aui = collections.defaultdict(list)
        self.atts = list()
        self.atts_by_code = collections.defaultdict(list)
        self.rank = list()
        self.rank_by_tty = collections.defaultdict(list)
        self.sty = list()
        self.sty_by_cui = collections.defaultdict(list)
        self.cui_roots = set()
        self.umls_version = umls_version
        self.store_atts = store_atts
        self.all_ont_codes = all_ont_codes
        self.ont_link = ont_link

    def load_tables(self):
        mrconso = UmlsTable.get("MRCONSO", self.data_root, 
                cacheable_sabs=self.all_ont_codes)
        mrconso_filt = {'SAB': self.ont_code, 'LAT': 'ENG'} 
        relev_cuis = set()
        for atom in mrconso.scan(filt=mrconso_filt):
            index = len(self.atoms)
            self.atoms_by_code[get_code(atom, self.load_on_cuis)].append(index)
            if not self.load_on_cuis:
                self.atoms_by_aui[atom[IDXS.MRCONSO.AUI]].append(index)
            self.atoms.append(atom)
            relev_cuis.add(atom[IDXS.MRCONSO.CUI])
        LOG.debug("length atoms: %d", len(self.atoms))
        LOG.debug("length atoms_by_aui: %d", len(self.atoms_by_aui))
        LOG.debug("atom example: %s", str(self.atoms[0]))

        mrconso_filt = {'SAB': 'SRC', 'CODE': 'V-%s' % self.ont_code} 
        for atom in mrconso.scan(filt=mrconso_filt):
            self.cui_roots.add(atom[IDXS.MRCONSO.CUI])
        LOG.debug("length cui_roots: %d" % len(self.cui_roots))


        mrrel = UmlsTable.get("MRREL", self.data_root,
                cacheable_sabs=self.all_ont_codes)
        gen_filt = {'SAB': self.ont_code} 
        field = IDXS.MRREL.AUI2 if not self.load_on_cuis else IDXS.MRREL.CUI2
        for rel in mrrel.scan(filt=gen_filt):
            index = len(self.rels)
            self.rels_by_aui_src[rel[field]].append(index)
            self.rels.append(rel)
        LOG.debug("length rels: %d", len(self.rels))

        mrdef = UmlsTable.get("MRDEF", self.data_root,
                cacheable_sabs=self.all_ont_codes)
        field = IDXS.MRDEF.AUI if not self.load_on_cuis else IDXS.MRDEF.CUI
        for defi in mrdef.scan(filt=gen_filt):
            index = len(self.defs)
            self.defs_by_aui[defi[field]].append(index)
            self.defs.append(defi)
        LOG.debug("length defs: %d", len(self.defs))

        if self.store_atts:
            mrsat = UmlsTable.get("MRSAT", self.data_root,
                    cacheable_sabs=self.all_ont_codes)
            field = IDXS.MRSAT.CODE if not self.load_on_cuis else IDXS.MRSAT.CUI
            for att in mrsat.scan(filt=gen_filt):
                index = len(self.atts)
                if not att[field]:
                    continue
                self.atts_by_code[att[field]].append(index)
                self.atts.append(att)
            LOG.debug("length atts: %d", len(self.atts))

        mrrank = UmlsTable.get("MRRANK", self.data_root,
                cacheable_sabs=self.all_ont_codes)
        for rank in mrrank.scan(filt=gen_filt):
            index = len(self.rank)
            self.rank_by_tty[rank[IDXS.MRRANK.TTY]].append(index)
            self.rank.append(rank)
        LOG.debug("length rank: %d", len(self.rank))

        mrsty = UmlsTable("MRSTY", self.data_root)
        for sty in mrsty.scan(filt=None):
            cui = sty[IDXS.MRSTY.CUI]
            if cui not in relev_cuis:
                continue
            index = len(self.sty)
            self.sty_by_cui[cui].append(index)
            self.sty.append(sty)
        LOG.debug("length sty: %d", len(self.sty))
        # Track the loaded status
        self.loaded = True
        LOG.info("%s tables loaded ...", self.ont_code)

    def terms(self):
        if not self.loaded:
            self.load_tables()
        # Note: most UMLS ontologies are 'load_on_codes' (only HL7 is load_on_cuis)
        for code in self.atoms_by_code:
            code_atoms = [self.atoms[row] for row in self.atoms_by_code[code]] 
            field = IDXS.MRCONSO.CUI if self.load_on_cuis else IDXS.MRCONSO.AUI 
            ids = map(lambda x: x[field], code_atoms)
            rels = list()
            for _id in ids:
                rels += [self.rels[x] for x in self.rels_by_aui_src[_id]]
            rels_to_class = list()
            is_root = False
            if self.load_on_cuis:
                rels_to_class = rels
                for rel in rels_to_class:
                    if rel[IDXS.MRREL.CUI1] in self.cui_roots:
                        is_root = True
                        break
            else:
                for rel in rels:
                    rel_with_codes = list(rel)
                    aui_source = rel[IDXS.MRREL.AUI2]
                    aui_target = rel[IDXS.MRREL.AUI1]
                    code_source = [ get_code(self.atoms[x], self.load_on_cuis) \
                                        for x in self.atoms_by_aui[aui_source] ]
                    code_target = [ get_code(self.atoms[x], self.load_on_cuis) \
                                        for x in self.atoms_by_aui[aui_target] ]
                    # TODO: Check use of CUI1 (target) or CUI2 (source) here:
                    if rel[IDXS.MRREL.CUI1] in self.cui_roots:
                        is_root = True
                    if len(code_source) != 1 or len(code_target) > 1:
                        raise AttributeError, "more than one or none codes"
                    if len(code_source) == 1 and len(code_target) == 1 and \
                        code_source[0] != code_target[0]:
                        code_source = code_source[0]
                        code_target = code_target[0]
                        # NOTE: the order of these append operations below is important.
                        # get_rel_code_source() - it uses [-1]
                        # get_rel_code_target() - it uses [-2]
                        # which are called from UmlsClass.toRDF().
                        rel_with_codes.append(code_target)
                        rel_with_codes.append(code_source)
                        rels_to_class.append(rel_with_codes)
            defs = [self.defs[x] for x in self.defs_by_aui[_id] for _id in ids]
            atts = [self.atts[x] for x in self.atts_by_code[code]]

            yield UmlsClass(self.ns, atoms=code_atoms, rels=rels_to_class,
                defs=defs, atts=atts, rank=self.rank, rank_by_tty=self.rank_by_tty,
                sty=self.sty, sty_by_cui=self.sty_by_cui,
                load_on_cuis=self.load_on_cuis, is_root=is_root, 
                ont_link=self.ont_link)

    def write_into(self, file_path, hierarchy=True):
        LOG.info("%s writing terms ... %s" % (self.ont_code, file_path))
        with open(file_path, "w", encoding=ENC) as fout:
            #nterms = len(self.atoms_by_code)
            fout.write(PREFIXES)
            comment = "RDF Version of the UMLS ontology %s; " +\
                      "converted with the UMLS2RDF tool " +\
                      "(https://github.com/ncbo/umls2rdf), "+\
                      "developed by the NCBO project."
            header_values = dict(
               label=self.ont_code,
               comment=comment % self.ont_code,
               versioninfo=self.umls_version,
               uri=self.ns
            )
            fout.write(ONTOLOGY_HEADER.substitute(header_values))
            for term in self.terms():
                fout.write(term.toRDF())

class UsageError(Exception):
    pass

def main():
    usage = "Usage: %prog [options] MEDLINE_RRF_DATA_DIR OUTPUT_DIR"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-u", "--umls-version", dest="umls_version",
            help="UMLS version to tag TTL files with", default="2012AB")
    parser.add_option("-s", "--source", dest="sources", action="append",
            help="Create TTL files for this source ontology. Can be "
            " specified multiple times; must be one of %s)" 
            % ', '.join(KNOWN_ONT_CODES))
    parser.add_option("--load-on-cuis", dest="load_on_cuis",
            default=False, action="store_true", 
            help="Load on CUIs (concept IDs) rather than "
            "AUIs (atom IDs) -- recommended for HL7 ")
    parser.add_option("-n", "--no-atts", dest="store_atts",
            default=True, action="store_false",
            help="Don't create triples for attributes (ie ignore MRSAT)")
    parser.add_option("--no-cache", dest="cache", action="store_false",
            default=True, help="Don't cache when creating multiple ontologies"
            " (saves memory but slows down, especially on smaller ontologies)")
    parser.add_option("--no-ont-link", dest="ont_link", action="store_false",
            default=True, help="Don't add in an explicit link to the source ontology "
            "using rdfs:isDefinedBy")
    (options, args) = parser.parse_args()
    try:
        try:
            data_root = args[0]
            output_dir = args[1]
        except IndexError:
            raise UsageError()
        if not path.isdir(output_dir) or not os.access(output_dir, os.W_OK):
            raise UsageError("Path %s does not exist or is not writable" % output_dir)
        for ps in options.sources:
            if ps not in KNOWN_ONT_CODES:
                raise UsageError("Unknown ontology %s" % ps)
    except UsageError:
        parser.print_help()
        raise

    output_file = os.path.join(output_dir, "umls_semantictypes.ttl")
    generate_semantic_types(data_root, STY_URL, output_file)
    if options.cache and len(options.sources) > 1:
        all_ont_codes = frozenset(options.sources)
    else:
        all_ont_codes = None
    for umls_code in options.sources:
        file_out = "%s.ttl" % umls_code
        output_file = path.join(output_dir, file_out)
        LOG.info("Generating %s (with load_on_cuis=%r, store_atts=%r, ont_link=%r)", umls_code, 
            options.load_on_cuis, options.store_atts, options.ont_link)
        ns = get_umls_url(umls_code)
        ont = UmlsOntology(umls_code, ns, data_root, 
                load_on_cuis=options.load_on_cuis,
                umls_version=options.umls_version,
                store_atts=options.store_atts,
                ont_link=options.ont_link,
                all_ont_codes=all_ont_codes)
        ont.load_tables()
        ont.write_into(output_file, hierarchy=(ont.ont_code != "MSH"))
        LOG.info("done!")
    LOG.info("generated MRDOC at global/UMLS level")


if __name__ == "__main__":
    main()
