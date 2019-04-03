# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 10:18:38 2017
"""
#import argparse
#import sys
from Bio.ExPASy import ScanProsite
import re
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


# Classe motif que contém a estrutura de dados
class motif:
    def __init__(self):
        self.score = 0
        self.pos = (0,0)   
        self.seq = ("")
        self.idSeq=("")
        self.acSeq=("")
        self.idMot=("")
        self.acMot=("")
        
    def setScore(self, score):
        self.score= score
        
    def setPos(self, start, stop):
        self.pos = (start, stop)
        
    def setSeq(self, seq):
        self.seq= seq
        
    def setId(self, idSeq):
        self.idSeq= idSeq
        
    def setAccesion(self, acSeq):
        self.acSeq= acSeq
        
    def setIdMot(self, idMot):
        self.idMot=idMot
        
    def setAcMot(self, acMot):
        self.acMot= acMot

#imprime no stdout o motif lido
"""parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfile", type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="output file to write the motifs")
args = parser.parse_args()"""

#Exemplos de input
"""
s="MGKNVVVLGTQWGDEGKGKIVDLLTQDAQVVVRYQGGHNAGHTLKINGVKTVLRLIPSGMLRPNVTCYIANGVVLSPQALLSEIKELEGNGINVRERLRISLACPLILPYHIALDKARETHMGKSAIGTTGRGIGPAYEDKVARRALRVGDLFHRDRFANKLTELLDYHNFVLTQYFKQPAVDLESLLGESLQWAEELRPMVCDVSACLHEHRKQGENILFEGAQGVYLDIDHGTYPYVTSSNTCVGSVINGAGFGPRYIDYVLGITKAYTTRVGGGPFPTELLDDVGKRIAERGQEFGAVTGRPRRCGWFDAVLLKRSIELNSISGLCVTKLDVLDGLEVLRIAVAYKDRDGNILSRPPLAADDFNDLLPVYEELPGWQESTADVTVMSDLPANARAYLKRIEEILGIPIDMLSTGPERDSTITLRGPFL"
s1="P98073"
s2= "ENTK_HUMAN"
"""

#Exemplo de motif
"""
mot= "[{'sequence_ac': 'USERSEQ1', 'start': 11, 'stop': 18, 'signature_ac': 'PS01266', 'level_tag': '(0)'}, {'sequence_ac': 'USERSEQ1', 'start': 133, 'stop': 144, 'signature_ac': 'PS00513', 'level_tag': '(0)'}]"
"""

#scanSeq() efetua o scan de todos os motifs existentes numa sequência, recorrendo a base de dados Prosite
def scanSeq(seq):
    keywords = {'CC':'/SKIP-FLAG=FALSE;'}
    #Pesquisa na base de dados ScanProsite
    response = ScanProsite.scan(seq, 'http://www.expasy.org', 'xml', **keywords )
    #Lê todos os motifs encontrados
    obj = ScanProsite.read(response)
    
    #Se necessário, imprime no terminal os motifs que encontrou
    #args.outfile.write(str(obj) + "\n")
    
    #Tendo já os motifs, irá passá-los como argumento quando chama a função parseMotif(), assim como a sequência de query
    res= (parseMotif(str(obj), seq))

    return res
   
#parseMotif() recebe um ou mais motifs contidos na sequência de query e preenche s estrutra com os dados corretos
def parseMotif(mot, seq):
    #se houver mais do que um motif, faz um split em "}, {" de maneira a dividir os motifs um por um
    motifs = re.split("}, {", mot)
    #faz o parser de cada item do motif. Por exemplo: "[{'sequence_ac': 'USERSEQ1'"
    myreg=re.compile(r"((\[(\{))?\'(\w+)\': ((\')?(\()?(\w+)(\()?(\')?)(\}\])?)")
    j=0
    res=[]
    #percorre todos os motifs
    for i in motifs:
        #Separa os "atributos" do motif pela virgula
        atributos= re.split(",", motifs[j])
        v=0
        ini=fim=0
        exp=motif()
        #percorre todos os atributos
        for m in atributos:
            #para cada atributo, faz o parse "myreg" e preenche na estrutura o valor de acordo com cada campo
            for data in re.findall(myreg, atributos[v]):
                if(data[3]=='score'): motif.setScore(exp,data[7]) 
                if(data[3]=='start'): ini=data[7]
                if(data[3]=='stop'): fim= data[7]
                if(data[3]=='sequence_id'): motif.setId(exp, data[7])
                if(data[3]=='sequence_ac'): motif.setAccesion(exp, data[7])
                if(data[3]=='signature_id'): motif.setIdMot(exp, data[7])
                if(data[3]=='signature_ac'): motif.setAcMot(exp, data[7])
            #no fim de percorrer todos os atributos, se ini e fim forem diferentes de zero, preenche o "campo" pos da estrutura    
            if(ini!=0 and fim!=0): motif.setPos(exp,ini,fim)
            i=int(ini)-1
            sequencia=""
            #Se a sequência tiver algum caracter que não seja maisculo ou se tiver números, ou seja,
            #se a sequência for em formato fasta e não o identificador ou acession, irá preencher o campo
            #seq com a sequência que corresponde ao motif
            if (seq.isalpha() and (not(any(char.isdigit() for char in seq)))):
                while(i<int(fim)):
                    sequencia= sequencia + seq[i]
                    i+=1
                motif.setSeq(exp, sequencia)
                
            v+=1
        #imprime todos os campos da estrutura
        print("Acession Sequência:", exp.acSeq)  
        print("Identificador Sequência:", exp.idSeq) 
        print("Acession Prosite Motif:", exp.acMot)  
        print("Identificador Prosite Motif:", exp.idMot)
        print("Posição:", exp.pos)
        print("Score:" ,exp.score) 
        if (seq.isalpha() and (not(any(char.isdigit() for char in seq)))):
            print("Sequência:", exp.seq)   
        res.append(exp)
        print("\n")
        j+=1
    
    return res
    
    


#Tentativa de fazer o scan na base de dados CDD
def scanCDD(seq):
    res= []
    response = NCBIWWW.qblast("blastp","cdd",seq,'https://blast.ncbi.nlm.nih.gov/Blast.cgi',format_type='XML')
    obj= NCBIXML.parse(response)
    
    for o in obj:
        alignments = sorted(o.alignments, key=lambda a: a.hsps[0].expect)[0:10]
        for a in alignments:
            res.append(a)
    return res


    