#! /usr/bin/env python

''' Phy-Bi
Written by Lachlan McGinness and Megan C McDonald
April 2015
MIT Licence

Program for collapsing phylogenetic trees into phylogenetic bi-trees.#
Phy-Bi is a program that allows biologists to gain a summary of the relatedness of different isolates.

Phy-Bi requires 3 inputs:
1.) A series of phylogenetic trees in Newick format. The filenames should take the form *.Genename
2.) A .txt file named 'genefiles.txt' which contains a the name of each file from 1 on a new line.
3.) A .txt file named 'isolatelist.txt' which contains the name of all isolates to be considered on the first line separated by commas.
## If a biologist only has the gene sequences for the different isolates, these can be converted into
phylogenetic trees by a RAxML (http://sco.h-its.org/exelixis/web/software/raxml/index.html) or a similar program.

Phy-Bi takes the phylogenetic tree corresponding to each gene and splits the isolates into two groups.
This two group pattern is called a bi-tree. For a given number of isolates there is a limited number
of possible bi-trees (combinations). For example if there are 16 possible combinations:
5 combinations which are 4-1
10 combinations which are 3-2
1 combination where all isolates are put into the same group. If all the isolates are  This is called a star. 

Phy-Bi then creates a file called 'BiTreesSummary.txt' which contains the frequency of occurance each bi-tree.
Plotting this data gives the biologist information about how closely related each isolate is.

#all - AllBin
One of the most useful functions of Phy-Bi is allowing biologists to find genes that match
a particular bitree pattern (which may correspond to experimental data). If the 'AllBin' function
is used then a file is created for each bi-tree pattern in a folder called 'genebins'.

#s- Search
Due to sequencing errors, occasionally genes can be missing from particular isolates. Phy-Bi will
ignore missing isolates and continue to operate, however this will lead to a larger number of possible
bi-tree combinations. Once the number of of isolates exceeds 6 it becomes impractical to find desired
bi-tree combinations by hand.
In these cases Phy-Bi has a search function. For Search to work the 'isolatelist.txt' input file must contain
a second line. This second line is a list of all isolates which the user wishes to be grouped together.
Phy-Bi will create a file named 'SearchResults.txt' which contains a list of each gene where every present
isolate in the second list is sorted into the same group.

#ab - Absence
The absence of certain genes may be of interest to the researcher. If the Absence function is used
then Phy-Bi will create a folder called 'AbsentBins', similar to genebins created by the AllBin
function. The files created in the AbsentBins folder will contain a list of all possible
combinations of absent genes. If both the Search and Absence functions are used then a second list
will be created in 'SearchResults.txt' which contains a list of all genes which have an absence
pattern which matches the second list inputted in the Search Function. It is recommented that
f is set to 0 if Absence is used.


There are a number of other minor functions

#d - Draw
Phy-Bi draws original phylogenic tree so user can compare output to the imput and check that the
program is giving reasonable results. Most useful to check on a small number of genes.

#f - Fraction
The user must give a number between 0 and 1. Phy-Bi will skip genes if too many isolates are missing.
This determines the minimum number of isolates that must be present for the gene not to be skipped.
If interested in Absence function, it is recommended that f is set to 0.

#i - Input directory
Specifies input directory. Default is current working directory.

#o - Output directory
Specifies output directory. Default is current working directory.

#v - Verbose
Gives more information about what PhyBi is doing. Can be used in combination with draw to check that Phy-Bi is giving the desired result.

#t - Threshold
Determines the distance requried for Phy-Bi to output a star. The default is 0.01 which corresponds
to a mutation rate of one in one hundred base pairs

#t2 - Threshold2
A number which indicates how closely related isolates need to be for Phy-Bi to automatically place
them in the same group. The default is 0.0002 which corresponds to a mutation rate of two in ten
thousand base pairs.
'''


### Packages
#import Bio
from Bio import SeqIO
from Bio import Phylo
import re
from operator import itemgetter
import os.path
import copy
import argparse
import shutil

### Argparses Start here
parser = argparse.ArgumentParser('TreeSorter3')
parser.add_argument('-ab', '--absence', action='store_true', help='If the Absence function is used then Phy-Bi will create a folder called AbsentBins, similar to genebins created by the AllBin function. The files created in the AbsentBins folder will contain a list of all possible combinations of absent genes. If both the Search and Absence functions are used then a second list will be created in SearchResults.txt which contains a list of all genes which have an absence pattern which matches the second list inputted in the Search Function. It is recommented that f is set to 0 if Absence is used.') # Change default to false
parser.add_argument('-all', '--allbin', action='store_true', help='This does what?') # Change default to false
parser.add_argument('-c', '--clear', action='store_true', help='Phy-Bi clears the "AbsenceBins" and "genebins" folders')
parser.add_argument('-d', '--draw', action='store_true', help='Phy-Bi draws original phylogenic tree so user can compare output to the imput and check that theprogram is giving reasonable results. Most useful to check on a small number of genes.')
parser.add_argument('-f', '--fraction', type=float, default=0.5, help ='A number between 0 and 1. This determines the minimum number of isolates that must be present for the gene not to be skipped. If interested in Absence function, it is recommended that -f is set to 0.')
parser.add_argument('-i', '--input', default=os.getcwd(), help='Path to input directory, [default: from local directory]')
parser.add_argument('-o', '--output', default=os.getcwd(), help='Path to place output, [default: to local directory]')
parser.add_argument('-v', '--verbose', action='store_true', help='Silences the default output comments given by Phy-Bi. Avoids cluttering the terminal.')
parser.add_argument('-s', '--search', action='store_true', help='For Search to work the isolatelist.txt input file must contain a second line. This second line is a list of all isolates which the user wishes to be grouped together. Phy-Bi will create a file named SearchResults.txt which contains a list of each gene where every present isolate in the second list is sorted into the same group.') # change default to false
parser.add_argument('-t', '--threshold', type=float, default=0.01, help ='requires a number which indicates how closely related isolates need to be for Phy-Bi to output a star, default is 0.01 which corresponds to a mutation rate of one in one hundred base pairs')
parser.add_argument('-t2', '--threshold2', type=float, default=0.0002, help ='requires a number which indicates how closely related isolates need to be for Phy-Bi to automatically place them in the same group, default is 0.0002 which corresponds to a mutation rate of two in ten thousand  base pairs')
args = parser.parse_args()

##### Inputs

# Isolates is a list that must be inputted by the user so Phy-Bi knows which isolates to consider

#genefiles.txt contains a list of all the files that Phy-Bi will act on. Each file name should be
#take the form *.genename (any text then a '.' then the name of the gene to be used by the biopthon
#package). Each gene should  seperated by a new line, either '\n' or '\r'
genestring = str(args.input + '/genefiles.txt')
filetoread = open(genestring, 'r')
genetext = filetoread.read()
genes = re.split('\n|\r', genetext)
genelist = []
for line in genes:
    if line != '':
        genelist.append(line)


#isolatelist.txt contains either one or two lines. There there are two lines they should be seperated
#by a new line, either '\n' or '\r'. Each line should contain a list of isolates seperated by a comma
        #and then a space ', '. The second line is only used if the Search function is used.
isolatestring = str(args.input + '/isolatelist.txt')
isolatefile = open(isolatestring, 'r')
isolatetext = isolatefile.read()
isolatelists = re.split('\n|\r', isolatetext)
if '' in isolatelists:
    isolatelists.remove('')
Isolatelist = re.split(', ', isolatelists[0])
Isolates = Isolatelist 
if len(isolatelists) >= 1:
    Searchlist  = re.split(', ', isolatelists[1])

### Outputs
if not os.path.exists(args.output):
    os.makedirs(args.output)
elif args.clear:
    if os.path.exists(os.path.join(args.output, 'AbsentBins')):
        shutil.rmtree(os.path.join(args.output, 'AbsentBins'))
    if os.path.exists(os.path.join(args.output, 'genebins')):
        shutil.rmtree(os.path.join(args.output, 'genebins'))
    

## The different Functions are contained below. These can be switched on or off depending on the
## User's Needs.
##### Functions


### AllBin
#AllBin = True    
comlist = []
comstring = []
if args.allbin:
    genebins_dn = os.path.join(args.output, 'genebins')
    if not os.path.exists(genebins_dn):
        os.makedirs(genebins_dn)

### Search
SearchTarget = []
if args.search:
    SearchTarget = Searchlist ### Will only find genes if it contains only these isolates
    # "St55", "St56", "St79", 'WAI326'
    SearchTarget.sort()
    SearchIterator = iter(SearchTarget)
    if args.verbose:
        print('Search for ')
        for target in SearchTarget:
            print(str(SearchIterator.next()))
    Matches = []
    Matches2 = []
    for isolate in SearchTarget:
        if isolate not in Isolates:
            print('Warning, one of the Isolates for which you are search is not included in the list of Isolates. It is likely that the search will return no results.')

        
### Absence
#Absence = True
if args.absence:
    abslist = []
    abstring = []
if args.absence:
    absencebins_dn = os.path.join(args.output, 'AbsentBins')
    if not os.path.exists(absencebins_dn):
        os.makedirs(absencebins_dn)



###MatrixMin is a function that will return the minimum value from a list of lists
def MatrixMin(A): 
    return (min(map(min, A)))

def DebugPrint(S):
    if args.verbose:
        print(S)

##### Outputs
#Output file list and output file names:
out_dn = args.output
searchpath = os.path.join(out_dn, 'SearchResults.txt')
writepath = os.path.join(out_dn, 'BiTreesSummary.txt') #Tree1[7:17] ### Writepath is the output file name. 



##### Phy-Bi
### This part of the program starts a loop so that Phy-Bi will run on every file in 'genefiles.txt'

if args.verbose:
    print('start')
    
writelist = []
AbSearch = []
for gene in genelist:
    gene2 = gene
    components = gene2.split('.')
    genename = components[1] ###Note that this is used later on 
    treenum = Phylo.read(os.path.join(args.input, gene), "newick")
    treenum.rooted = True
    ##If you wish to see what the raw phylogenic trees look like you can use the code below:

    if args.draw:
        Phylo.draw_ascii(treenum)
    ##If you look at a small set of trees and the program output, you can check that you agree with decisions that Phy-Bi makes.

    ### Define the names of the isolates (referred to as 'islets' for short).
    n = 1
    islets = []
    islets1 = []
    for Isolate in Isolates:
        islets1.append([str(Isolate + "-" + genename), n, Isolate])
        n = n + 1

    #If isolates on the list to be searched are absent for a particular gene, this block of code removes them temporarily
    #if args.search:
    SearchTargetTemp = copy.deepcopy(SearchTarget)
    AbsentGenes = []
    for islet in islets1:
        x = 1
        try: 
            treenum.find_clades(name = islet[0]).next()
            islets.append(islet)
        except StopIteration:
            if islet[2] in SearchTargetTemp:
                        SearchTargetTemp.remove(islet[2])
            AbsentGenes.append(islet[2])
            x = 0

    
    if args.absence: 
       if args.search:
            AbsentGenes.sort
            Abtemplist = []
            if SearchTargetTemp == AbsentGenes: #(note less than (<=) gives you the subset)
                for islet in islets:
                    Abtemplist.append(islet[2])
                AbSearch.append([genename,  ''.join(AbsentGenes), ''.join(Abtemplist)])

                
    
    ### This loop is designed to check if all the islets are so closely related that we will 
    ### just throw them together as one group. If they are all within a certain cut of distance
    ### then the program will just output a 'star'. Given that most genes are so similar that 
    ### most will be stars it was important to do this check early to save time.

    ###     
    if args.search:
        if args.verbose:
            print(len(islets), len(islets1))

    ### The code below determines how many isolates can be missing from a gene before Phy-Bi skips the gene.
    ### -f command
    #if True:   #Run regardless of how many isolates are present
    #if  len(islets) == len(islets1):  ###Run only if all isolates are present
    if len(islets) >= len(islets1)*args.fraction:  ###Run if at least half the isolates are present
        distances = []
        distances2 = []
        n = 1
        while n <= len(islets):
            m = n
            while m <= len(islets):
                distances.append(treenum.distance(islets[n-1][0], islets[m-1][0]))
                distances2.append([str(str(n)+'_'+str(m)), treenum.distance(islets[n-1][0], islets[m-1][0])])
                m = m+1
            n = n+1
        
        writelisttemp = []
        ### Arbitary cut off point to determine closeness of relation
        ### -t command
        if max(distances) < args.threshold: # default is one basepair in 100
            writelisttemp.append('star')
            if args.verbose:
                print (writelisttemp[0])
            writelist.append('star') 
            writelist.append('\n')
            if True:
                if 'star' in comstring:
                    for item in comlist:
                        if item[0] == 'star':
                            item.append(genename)
                else:
                    comstring.append('star')
                    comlist.append(['star', genename])
        else: 
            ###### Part 1 ######
            ### This loop creates groups of islets that are very closely related. By doing this here it 
            ### will reduce the number of steps in part 2 and part 3 which are more computationally 
            ### expensive. Overall this should decrease total running time.
            groups = [] ### Will be a list of all the groups
            n = 1
            for islet in islets:
                i = 0
                m=0 #A vairable used to check whether islet1 is already used
                name = str('group'+str(n)) # For example, group1, group2 e.t.c
                for group in groups:
                    if islet in group: # This tells the program if the islet is already in a group.
                        i = 1 
                if i == 0: # If the islet is already in a group then the next code won't run
                    for group in groups: # If the islet isn't in a group then it will check to see if
                        # it is closely related to one of the other groups.
                        if treenum.distance(group[1][0], islet[0]) <= args.threshold2: #arbirary cut off for very closely related individuals
                            group.append(islet) # Adds the islet to the closely related group
                            m = 1
                    if m == 0: #If the islet was not closely related to any of the other groups it starts a new group.
                        groups.append([name, islet])
                        n = n + 1

            ####### Part 2 #####
            ### This loop tells each group how closely related it is to each other group. This information
            ### is stored at the end of the group. Each group takes the form:
            ### [member1, member2, ... [X, x1, x2...]] Where: X is a number indicating on average how 
            ### closely related you are to all the other groups, and x1, x2, ... are numbers indicating how
            ### closely related this group is to group1, group2, ...
            xsTab = []
            for group1 in groups: 
                xs = [] ### xs is a list of variables which measures the distance (x) between each group
                for group2 in groups:
                    if group1 != group2:
                        xs.append(treenum.distance(group1[1][0], group2[1][0]))
                    else: 
                        xs.append(1)
                group1.append(xs)
                xsTab.append(xs)

            #####  Part 3 ####
            ### This loop merges all the groups until only two groups remain
            ### It starts out by finding the group which is most closely related to 
            ### all the other groups. (3a)
            ### It then identifies the most closely related group. (3b)
            ### It then moves the members into the new group. (3c)
            groups2 = copy.deepcopy(groups)
            numberofgroups = len(groups2)
            while numberofgroups >= 3:
                numberofgroups = len(groups2)
                groups3 = copy.deepcopy(groups2)
                b = MatrixMin(xsTab)
                for row in xsTab:
                    if b in row:
                        IndexR = row.index(b)
                        IndexA = xsTab.index(row)
                        break
                length = len(groups3[IndexR])-1
                toinsert = groups3[IndexR][1:length] #This ensures that you move the isolates but not 
                #the group name at the beginning or the numbers at the end
                for group in groups2:
                    if groups3.index(group) == IndexA: #(group that it is being appended to)
                        for item in toinsert:
                            group.insert(1, item)
                        for entry in group[-1]: ###Deals with all other numbers
                            group[-1][group[-1].index(entry)] = (groups3[IndexR][-1][group[-1].index(entry)]*(len(groups3[IndexR])-2) + entry*(len(groups3[IndexA])-2))/(len(groups3[IndexR])-2 + len(groups3[IndexA])-2) ### THIS Needs to be weighted according to number of isloates
                        group[-1][IndexA]= 1 # fixes the self co-ordinate
                        group[-1].remove(group[-1][IndexR]) # removes group which is to be removed  
                    elif groups3.index(group) != IndexR: 
                        if groups3.index(group) != IndexA: #(All other groups)
                            group[-1][IndexA] = (group[-1][IndexA]*(len(groups3[IndexA])-2) + group[-1][IndexR]*(len(groups3[IndexR])-2))/((len(groups3[IndexA])-2) + (len(groups3[IndexR])-2))### Needs to be weighted according to number of isolates
                            group[-1].remove(group[-1][IndexR]) 
                groups2.remove(groups2[IndexR])
                numberofgroups = numberofgroups -1
                numberofgroups = len(groups2)
                xsTab = []
                for group in groups2:
                    xsTab.append(group[-1])
                
            
            ####### Part 4 #####
            
            ###This section simply converts to Newick Format
            for group in groups2: 
                del group[0] ### Removes 'group1' at beginning of group
                del group[-1] ### Removes all the numbers at the end of the group
                
            ### This section of code gets it to sort the islets so they always output in a consistent order
            newgroups = []
            templist1 = [(islet[1], islet) for islet in groups2[0]]
            templist2 = [(islet[1], islet) for islet in groups2[1]]
            templist1.sort()
            templist2.sort()
            newgroup1 = [group[1] for group in templist1]
            newgroup2 = [group[1] for group in templist2]
            ### This section of code ensures that the largest group is always listed first.
            if len(newgroup1) >= len(newgroup2):
                newgroups = [] 
                newgroups.append(newgroup1)
                newgroups.append(newgroup2)
            else:
                newgroups = [] 
                newgroups.append(newgroup2)
                newgroups.append(newgroup1)
                
        
            
            writelisttemp = []
            ### writing to output file in Newick Format
            n = 0 #n Counts the groups (used to ensure there is no comma for the last group)
            writelist.append('(')
            writelisttemp.append('(')
            for newgroup in newgroups:
                writelist.append('(')
                writelisttemp.append('(')
                m = 1 #m counts the items in the groups (to ensure there is no comma for the last item)
                for item in newgroup:
                    writelist.append(item[2])
                    writelisttemp.append(item[2])
                    if m != len(newgroup):
                        writelist.append(',')
                        writelisttemp.append(',')
                    m = m +1
                writelist.append(')')
                writelisttemp.append(')')
                if n <= 0:
                    writelist.append(',')
                    writelisttemp.append(',')
                    n = n + 1
                else:
                    writelist.append(')')
                    writelisttemp.append(')')
            writelist.append('\n')
            
            treename = ''.join(writelisttemp)
            if True:
                if treename in comstring:
                    comindex = comstring.index(treename)
                    comlist[comindex].append(genename)
                else:
                    comstring.append(treename)
                    comlist.append([treename, genename])    

            # This prints the classificiation so the user can see Phy-Bi's decision for each gene.
            sent_str = ""
            for i in writelisttemp:
                sent_str += str(i) 
            sent_str = sent_str[:-1]
            if args.verbose:
                print(sent_str)
            

            #### This part of the program takes AbsentGenes and uses it to create a list an absence pattern if this pattern is not already on ablist. It appends the gene to that absence pattern.
            if args.absence:
                AbsentGenes.sort()
                AbsencePattern = ''.join(AbsentGenes)
                if AbsencePattern in abstring:
                    absindex = abstring.index(AbsencePattern)
                    abslist[absindex].append(genename)
                else:
                    abstring.append(AbsencePattern)
                    abslist.append([AbsencePattern, genename])
                


            if args.search:
                templist = []
                templist2 = []
                for item in newgroups[1]:
                    templist.append(item[2])
                templist.sort()
                for item in newgroups[0]:
                    templist2.append(item[2])
                templist.sort()
                templist2.sort()
                if templist == SearchTargetTemp:
                    Matches.append(genename)
                    Matches2.append(genename + '\t' + str(len(SearchTargetTemp)) + '\t' + str(len(islets))+ '\t' + ''.join(templist) +'\t' + ''.join(templist2)) 
                if templist2 == SearchTargetTemp:
                    Matches.append(genename)
                    Matches2.append(genename + '\t' + str(len(SearchTargetTemp)) + '\t' + str(len(islets))+ '\t' + ''.join(templist2) + '\t' + ''.join(templist))

if args.verbose:
    print('Writing every found combination.')

for item in comlist:
    with open(writepath, 'w') as out:
        for item in comlist:
            out.write(item[0])
            out.write('\t')
            out.write(str(len(item)-1))
            out.write('\n')


### Creates a list of all the absence patterns
if args.absence:
    if args.verbose:
        print('Writing Absence Patterns')
    for combination in abslist:
        if combination[0] != '':
            outpath = os.path.join(absencebins_dn, combination[0] + '.txt')
        #outpath = str('AbsentBins/Absent' + combination[0] +'.txt') #
            with open(outpath, 'w') as out:
                for item in combination:
                    out.write(item)
                    out.write('\n')

if args.search:
    if args.verbose:
        print('Writing Matching Gene Combinations')
    with open(searchpath, 'w') as out:
        out.write('List of Genes with desired grouping')    
        out.write('\n') 
        out.write('Gene' + '\t' + '# Desired Isolates' + '\t' + '# Isolates Present ' + '\t' + 'Group1' + '\t' + 'Group2')
        out.write('\n') 
        for match in Matches2:
            out.write(match)
            out.write('\n')
        if args.absence:
            out.write('\n')
            out.write('\n')
            out.write('\n')
            out.write('\n')
            out.write('\n')
            out.write('Matching Absence Patterns')
            out.write('\n')
            out.write('Gene' + 'Absent Isolates' + 'Present Isolates')
            out.write('\n')
            for line in AbSearch:
                for item in line:
                    out.write(item)
                    out.write('\t')
                out.write('\n')

     
### Creates lists of genes for the trees
if args.allbin: 
    if args.verbose:
        print('Writing all bins.') 
    for combination in comlist:
        binfilename = combination[0].replace('),(', '_')
        binfilename = binfilename.replace(',', '')
        binfilename = binfilename.replace('(', '')
        binfilename = binfilename.replace(')', '')
        outpath = os.path.join(genebins_dn, binfilename + '_genes.txt')
        with open(outpath, 'w') as out:
            for item in combination:
                out.write(item)
                out.write('\n')

if args.verbose:
    print('done')
