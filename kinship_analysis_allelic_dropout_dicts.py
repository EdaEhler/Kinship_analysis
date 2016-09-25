###############################
# Kinship analysis using SNPs #
# --------------------------- #
# Edvard Ehler, Ph.D.,        #
# Institute of Anthropology,  #
# UAM Poznan, Poland,         #
# 2016                        #
# eda.ehler@seznam.cz         #
###############################

# Based on Fung&Hu (2008) Statistical DNA Forensics: Theory, Methods and Computation
"""
For two persons X and Y, define the relatedness coefficients (k 0 ,2k 1 ,k 2 ) as
k0 = P (neither allele of X is identical by descent to alleles of Y);
k1 = P (one or the other of the alleles of X is ibd to one of the alleles of Y, but the second allele is not); 
k2 = P (both alleles of X are ibd to those of Y).
"""
#--------------------
# IMPORTS
#--------

from collections import namedtuple  # používám k volání příbuzenských koeficientů
from math import pow                # cca 15% rychlejší než klasickej pow
from random import random


# VARIABLES
#----------

snp_list = []

# {'rs112' : {"G":0.3, "T":0, "C":0.7, "A":0}, 'rs2341': {"G":0.2, "T":0.8, "C":0, "A":0}}
freq_dict = {}


# table 3.13, chapter 3.6, page 43
relationship = namedtuple("rel_coefficients", "k0, k1, k2")

#----
parentChild = relationship(0, 1, 0)
fullSiblings = relationship(0.25, 0.5, 0.25)
halfSiblings = relationship(0.5, 0.5, 0)
grandparentGrandchild = relationship(0.5, 0.5, 0)
uncleNephew = relationship(0.5, 0.5, 0)
firstCousins = relationship(0.75, 0.25, 0)
secondCousins = relationship(0.9375, 0.0625, 0)
unrelated = relationship(1, 0, 0)



SNP_INFO = "SNP_allele_freqs.csv"
SAMPLES_GENOTYPES = "3samples_genotypes.csv"
# dle Anny - P(false positive) of homozygotes
# všechny homozygoty ve funkci divide_et_impera budu testovat, zda nejsou false positive,
# jestli jo, tak je budu brát jako heterozygoty
ALLELIC_DROPOUT = 0.00159
ALLELIC_DROPOUT_PROBS = "3samples_allelic_dropouts.csv"



# na vzorky
kz1 = {}
kz5 = {}
kz4 = {}

# na P(false homozygote) slovník (= allelic drop-out = ado)
kz1_ado = {}
kz5_ado = {}
kz4_ado = {}

snp_counter_nonzero = 0
snp_counter = 0


#-----------------
# LOADING SNP info + allel frequencies + samples genotypes
#-------------------------------------

with open(SNP_INFO, mode="r", encoding="utf-8") as snpIN:
    # načti do dvou slovníků, v jednom budou pouze názvy rsxXX(asi tuple), druhý slovník bude odkazovat na jejich parametry
    for radek in snpIN:
        radek = radek.strip().split(";")
        
        # jméno do snp_listu, abych po něm mohl pak cyklit
        snp_list.append(radek[0])
        
        # frekvence alel do freq_dict
        # nejdříve však defaultní hodnoty
        freq_dict[radek[0]] = {"G":0, "T":0, "C":0, "A":0}
        freq_dict[radek[0]][radek[1]] = radek[3]
        freq_dict[radek[0]][radek[2]] = radek[4]

with open(SAMPLES_GENOTYPES, mode="r", encoding="utf-8") as genoIN:
    # načtu genotypy vzorků
    for radek in genoIN:
        # nechci vykomentované řádky a N genotypy
        if not radek.startswith("#"):
            if not "N" in radek:
                radek = radek.strip().split(";")
                #h2[radek[0]] = radek[1]
                #h4[radek[0]] = radek[2]
                kz1[radek[0]] = radek[1]
                kz5[radek[0]] = radek[2]
                kz4[radek[0]] = radek[3]

with open(ALLELIC_DROPOUT_PROBS, mode="r", encoding="utf-8") as adIN:
    # načtu genotypy vzorků
    for radek in adIN:
        # nechci vykomentované řádky a N genotypy
        if not radek.startswith("#"):
            if not "N" in radek:
                radek = radek.strip().split(";")
                #h2[radek[0]] = radek[1]
                #h4[radek[0]] = radek[2]
                kz1_ado[radek[0]] = float(radek[1])
                kz5_ado[radek[0]] = float(radek[2])
                kz4_ado[radek[0]] = float(radek[3])


        
# FUNCTIONS
#----------
def divide_et_impera(snp, genotype1, genotype2, alleleCount, scenario=parentChild, jmeno1="prvni", jmeno2="druhy"):
    # dle poměru alel v genotypech zavolá odpovídající funkci
    # snp - name of SNP, string ('rs12345')
    # genotype1 - first individual genotype, string ("CC")
    # genotype2 - second individual genotype, string ("AC")
    # alleleCount - number of alleles across all loci in the 2 individuals tested, int (124)
    # scenario - relationship namedtuple defined earlier with relatedness coefficients (k0, k1(which is in fact 2k1, but naming problems made me to name it just k1), k2)
    # jmeno1, jmeno2 - jmeno vzorku, dle toho zařídím slovník pro allelic dropout    
    global snp_counter, snp_counter_nonzero
    
    #------------
    # blok definice allelic drop-out slovníku (ado)
    # dle toho, co přijde za jméno do funkce, volím slovník
    ado1 = {}
    ado2 = {}
    if jmeno1.upper() == "KZ1":
        ado1 = kz1_ado
    elif jmeno1.upper() == "KZ4":
        ado1 = kz4_ado
    else:
        print("jmeno1:", jmeno1)
        raise NameError("jmeno1 has unknown value (not KZ1, KZ4).")
        
    if jmeno2.upper() == "KZ4":
        ado2 = kz4_ado
    elif jmeno2.upper() == "KZ5":
        ado2 = kz5_ado
    else:
        print("jmeno2:", jmeno2)
        raise NameError("jmeno2 has unknown value (not KZ4, KZ5).")
    #-------------------
    
    
    # pomocná proměnná na testování rozřazovacího algoritmu
    branch = ""
    
    #-------------------
    #Rozřazování dle genotypů:
    # AA, AA
    if (genotype1 == genotype2) and (genotype1[0] == genotype1[1]):
        branch = "aaaa"
        allele1 = genotype1[0]
        allele2 = genotype1[0]
        # rozstřel genotypů na allelic dropout
        drop_out_roll = random()
        if drop_out_roll <= ado1[snp] * ado2[snp]: # pravděpodobnost, že jsou oba false positive
            funkce = ab_ab
        elif drop_out_roll <= ado1[snp]: # pravděpodobnost, že je jeden false positive
            funkce = aa_ab
        elif drop_out_roll <= ado2[snp]: # pravděpodobnost, že je druhý false positive
            funkce = aa_ab
        else:
            funkce = aa_aa
    
    # AB, AB
    elif (genotype1 == genotype2) and (genotype1[0] != genotype1[1]):
        branch = "abab"
        allele1 = genotype1[0]
        allele2 = genotype1[1]
        # rozstřel genotypů na allelic dropout
        drop_out_roll = random()
        # první možnost nedám - to by znamenalo, že se oba mohou změnit oběma směrama
        #if drop_out_roll <= ado1[snp] * ado2[snp]: # pravděpodobnost, že jsou oba false positive
        #    funkce = ab_ab
        if drop_out_roll <= ado1[snp]: # pravděpodobnost, že je jeden false positive
            funkce = aa_ab
        elif drop_out_roll <= ado2[snp]: # pravděpodobnost, že je druhý false positive
            funkce = aa_ab
        else:
            funkce = ab_ab
    
    # AA, BB
    elif (genotype1 != genotype2) and (genotype1[0] == genotype1[1]) and (genotype2[0] == genotype2[1]):
        branch = "aabb"
        allele1 = genotype1[0]
        allele2 = genotype2[0]
        # rozstřel genotypů na allelic dropout
        drop_out_roll = random()
        if drop_out_roll <= ado1[snp] * ado2[snp]: # pravděpodobnost, že jsou oba false positive
            funkce = ab_ab
        elif drop_out_roll <= ado1[snp]: # pravděpodobnost, že je jeden false positive
            funkce = aa_ab
        elif drop_out_roll <= ado2[snp]: # pravděpodobnost, že je druhý false positive
            funkce = aa_ab
        else:
            funkce = aa_bb
    
    # AA, AB
    elif (genotype1 != genotype2) and (genotype1[0] == genotype1[1]) and (genotype2[0] != genotype2[1]):
        branch = "aaab"
        allele1 = genotype1[0]
        # nevím, jestli mi přijde genotype2 AB nebo BA
        allele2 = genotype2[1] if genotype2[1] != genotype1[0] else genotype2[0]
        # rozstřel genotypů na allelic dropout
        drop_out_roll = random()
        #if drop_out_roll <= ado1[snp] * ado2[snp]: # pravděpodobnost, že jsou oba false positive
        #    funkce = ab_ab
        if drop_out_roll <= ado1[snp]: # pravděpodobnost, že je jeden false positive
            funkce = ab_ab
        #elif drop_out_roll <= ado2[snp]: # pravděpodobnost, že je druhý false positive
        #   funkce = aa_bb
        else:
            funkce = aa_ab
    
    # AB, AA
    elif (genotype1 != genotype2) and (genotype1[0] != genotype1[1]) and (genotype2[0] == genotype2[1]):
        branch = "abaa"
        allele1 = genotype2[0]
        # nevím, jestli mi přijde genotype1 AB nebo BA
        allele2 = genotype1[1] if genotype1[1] != genotype2[0] else genotype1[0]
        # rozstřel genotypů na allelic dropout
        drop_out_roll = random()
        #if drop_out_roll <= ado1[snp] * ado2[snp]: # pravděpodobnost, že jsou oba false positive
        #    funkce = ab_ab
        #elif drop_out_roll <= ado1[snp]: # pravděpodobnost, že je jeden false positive
        #    funkce = aa_ab
        if drop_out_roll <= ado2[snp]: # pravděpodobnost, že je druhý false positive
            funkce = ab_ab
        else:
            funkce = aa_ab
            
    # frekvence alel ve srovnávací populaci (bráno z ensemblu GRCh37)
    f1 = float(freq_dict[snp][allele1])
    f2 = float(freq_dict[snp][allele2])
    
    # test prints - byly špatné indexy v if-else bloku - už jsou OK
    """
    print(branch)
    print("genotyp:", genotype1, genotype2)
    print("allele:", allele1, allele2)
    print("pi,pj:", f1, f2)
    print("P(ano):", funkce(f1, f2, koef=scenario))
    print("P(ne):", funkce(f1, f2, koef=unrelated))
    print("LR:", funkce(f1, f2, koef=scenario) / funkce(f1, f2, koef=unrelated))
    input()
    """
    
    likelihoodRatio = funkce(f1, f2, koef=scenario) / funkce(f1, f2, koef=unrelated)
        
    snp_counter += 1
    
    if likelihoodRatio == 0:
        #print('zero', snp)
        # děje se zejména při parent-child scénáři, když se neshodují genotypy
        # pravděpodobnost mutací nebo silent alleles (Pinto et al. 2013, FSI:Genetics) -> 0.001-0.005
        #----------------
        # dle Borsting et al. 2011 (FSI:Genetics) počítá u rozdílných homozygotů jaby by se tam
        # objevila "silent allele" (třeba nějaká technická chyba, že ji nenašli).
        # pravděpodobnost silent allele je 1/(n+1), kde n = počet alel na všech lokusech u těchto dvou individuí
        # vynásobeno konzervativním odhadem pravděpodobnosti mutace u SNPů = 10E-6
        
        # print("+++zero+++")
        #return 0.000001 * (1/(alleleCount + 1))
        return 0
        #print("uvnitr divide_et_impera:", jmeno1)
        
        
    else:
        snp_counter_nonzero += 1
        return likelihoodRatio
    

# Fung&Hu 2008, table 5.1, page 80
# funkce, které počítají pravděpodobnost joint genotype probability za předpokladu HWE
# vstupují do nich frekvence alely 1 a 2 (f1, f2) a příbuzenské koeficienty, dle použitého scénáře
# výstup je P(Z|Y, H) - pravděpodobnost, že genotypy Z a Y mají alely identical-by-descend (ibd),
# za předpokladu hypotézy (scénáře) H (třeba že jsou siblings, nebo unrelated, nebo uncle-nephew...)
def aa_aa(f1, f2, koef):
    return koef.k0 * pow(f1, 4) + koef.k1 * pow(f1, 3) + koef.k2 * pow(f1, 2)

def aa_ab(f1, f2, koef):
    return 2 * koef.k0 * pow(f1, 3) * f2 + koef.k1 * pow(f1, 2) * f2


def aa_bb(f1, f2, koef):
    return koef.k0 * pow(f1, 2) * pow(f2, 2)
    
def ab_ab(f1, f2, koef):
    #print("abab:", (4 * koef.k0 * pow(f1, 2) * pow(f2, 2)) + (koef.k1 * pow(f1,2) * f2) + (koef.k1 * f1 * pow(f2, 2)) + (2 * koef.k2 * f1 * f2))
    return (4 * koef.k0 * pow(f1, 2) * pow(f2, 2)) + (koef.k1 * pow(f1,2) * f2) + (koef.k1 * f1 * pow(f2, 2)) + (2 * koef.k2 * f1 * f2)
    
#---------------------

def allele_count(sample1, sample2):
    
    allele_n = 0
    
    for i in sample1:
        geno1 = sample1[i]
        geno2 = sample2[i]
        # zanedbávám 3 a více alel, pouze bi-alelické lokusy
        allele_n += 1 if geno1 == geno2 and geno1[0] == geno1[1] else 2
        #print(i, allele_n)
    
    return allele_n
    
    
    
#---------------------
# SKRIPT
#-------

# kz4 vs kz5



def run_kinship_analysis(sample1, sample2, hypothesis, hypothesis_name, alleleCount, name1='prvni', name2='druhy', repeats=100):
    # přidán parametr name1,name2 - jméno pro výběr správného allelic drop-out slovníku
    
    # pomocná funkce na projetí všech definovaných kombinací vstupních parametrů
    # přidělal jsem možnost opakování výpočtu pro případ silent allele, allelic drop-in/drop-out
    # idea je taková, že provedu výpočet 1000x-100 000x a vezmu průměr
    # možná by byl lepší resampling??
    
    global snp_counter, snp_counter_nonzero
    
    snp_counter = 0
    snp_counter_nonzero = 0
    result = 1
    
    result_list = []
    
    # Opakování
    #--------------
    for _ in range(repeats):
        result = 1      # při každé rundě si vynuluju vysledek
        for i in sample1:
            
            #print(i, 'result = ', result)
            try:
                result *= divide_et_impera(i, sample1[i], sample2[i], alleleCount, scenario=hypothesis, jmeno1=name1, jmeno2=name2)
            except IndexError:
                #print(i, kz1[i])
                #print(i, kz5[i])
                pass
            #print("--")
        result_list.append(result)
    #--------------
    
    # zprůměrování výsledku
    result = sum(result_list)/repeats
    
    print(name1 + " vs. " + name2)
    print("Scenario:", hypothesis_name + ",", hypothesis)    
    print("Likelihood Ratio (p(scenario)/p(unrelated)):", result)
    print("Bayes. estimate of probability of the scenario (prior probability = 0.5):", str(round((result/(result + 1))*100, 5)) + "%")
    print("SNPs tried:", snp_counter/repeats)    
    print("SNPs with non-zero result:", snp_counter_nonzero/repeats)
    print("----------------------------------------------------")
    print()
    #input()


scenarios_bag = (parentChild, fullSiblings, halfSiblings, grandparentGrandchild, uncleNephew, firstCousins, secondCousins)
scenarios_names = ('parent-child', 'full-siblings', 'half-siblings', 'grandparent-grandchild', 'uncle-nephew', 'first cousins', 'second cousins')

#----------

pocet_alel = allele_count(kz4, kz5)
print("Allele count:", pocet_alel)
    
for n, hypo in enumerate(scenarios_bag):
    run_kinship_analysis(kz4, kz5, hypo, scenarios_names[n], pocet_alel, name1='KZ4', name2='KZ5', repeats=10000)
    #input()



print("============================================")
print("Algorithm loops count: 10000")
print("Allelic drop-out check - using dictionary of P(false allele) unique for each SNP for each sample.")
print("No silent-allele correction, just return 0 in case of opposite homozygotes with no drop-out.")
print("********************************************")
  

    
    


    
