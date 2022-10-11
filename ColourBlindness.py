
import numpy as np
import scipy as sp # may be useful to compute probabilities
import time # may be useful to check the execution time of some function
from tqdm import tqdm

"""
Please refer to lecture slides.
Please refer to README file.
All the functions that you define must be able to handle corner cases/exceptions
"""

"""
Starting and ending locations (indices) of red and green exons in the reference sequence - Begins

1. Red Exon Locations
"""
RedExonPos = np.array([
    [149249757, 149249868], # R1
    [149256127, 149256423], # R2
    [149258412, 149258580], # R3
    [149260048, 149260213], # R4
    [149261768, 149262007], # R5
    [149264290, 149264400]  # R6
    ])
"""
2. Green Exon Locations
"""
GreenExonPos = np.array([
    [149288166, 149288277], # G1
    [149293258, 149293554], # G2
    [149295542, 149295710], # G3
    [149297178, 149297343], # G4
    [149298898, 149299137], # G5
    [149301420, 149301530]  # G6
    ])
"""
Starting and ending locations (indices) of red and green exons in the reference sequence - Ends
"""    
def loadLastCol(filename):

    # function body - Begins

    file = open(filename,'r')    # This file extst
    # LastCol = file.readlines()
    # LastCol=[x[:-1] for x in LastCol]
    # LastCol="".join(LastCol)
    loadLastCol = file.read()
    file.close()


    # function body - Ends
    return loadLastCol #string data type

def loadRefSeq(filename):
    """
    Input: Path of the file containing the reference sequence.
    Output: The reference sequence in string format.
    """
    # function body - Begins

    #Read reference
    file = open(filename,'r')      # This file extst
    # file.readline()
    # RefSeq = file.readlines()
    # RefSeq=[x[:-1] for x in RefSeq]
    # RefSeq="".join(RefSeq)
    RefSeq = file.read()

    file.close()

    

    # function body - Ends
    return RefSeq # string data type

def loadReads(filename):
    """
    Input: Path of the file containing all the reads.
    Output: A list containing all the reads.
    """
    # function body - Begins

    file = open(filename,'r')
    file.readline()
    Reads = file.readlines()
    file.close()

    # function body - Ends
    return Reads # list of strings

def loadMapToRefSeq(filename):
    """
    Input: Path of the file containing mapping of the first column to the reference sequence.
    Output: numpy integer array containing the mapping.
    """
    # function body - Begins

    #Read indices
    file = open(filename,'r')     # This file extst
    MapToRefSeq = file.readlines()
    MapToRefSeq=np.array(MapToRefSeq,dtype=int)
    file.close()

    # function body - Ends
    return MapToRefSeq # numpy integer array


def MatchReadToLoc(read):
    """
    Input: a read (string)
    Output: list of potential locations at which the read may match the reference sequence. 
    Refer to example in the README file.
    IMPORTANT: This function must also handle the following:
        1. cases where the read does not match with the reference sequence
        2. any other special case that you may encounter
    """

    def get_Rel_address(add):
        main=add//100
        rel=add%100
        return main,rel

    def get_loc(read):

        def get_Bias(char):
            if char=='A':
                Bias =0
            elif char=='C':
                Bias = num_of_A
            elif char=='G':
                Bias = num_of_A+num_of_C
            elif char=='T':
                Bias = num_of_A+num_of_C+num_of_G
            return Bias
     
        #Initial Condition
        lower_main_adderess=0
        upper_main_adderess=num_of_A+num_of_C+num_of_G+num_of_T
        upper_main_adderess+=1
        read_len=len(read)

        #Searching in reverse for each suffix
        for i in reversed(range(read_len)):
            Char = read[i]
            upper_1= Milestone[Char][get_Rel_address(lower_main_adderess)[0]]
            lower_1=bwt[get_Rel_address(lower_main_adderess)[0]][get_Rel_address(lower_main_adderess)[1]:].count(Char)-1
            #Searching for the last C 
            upper_2=Milestone[Char][get_Rel_address(upper_main_adderess)[0]]
            lower_2=bwt[get_Rel_address(upper_main_adderess)[0]].count(Char)-bwt[get_Rel_address(upper_main_adderess)[0]][:get_Rel_address(upper_main_adderess)[1]+1].count(Char)

            if (upper_1-lower_1)-(upper_2-lower_2)>0:
                return lower_main_adderess,upper_main_adderess,i+1
            else:
                lower_main_adderess=upper_1-lower_1-1+get_Bias(Char)
                upper_main_adderess=upper_2-lower_2-1+get_Bias(Char)
        return lower_main_adderess,upper_main_adderess,0

    List_of_matches=[]
    read = read[:-1].replace('N','A')
    rev=read[::-1]
    replacements=[('A','T'),('G','C'),('T','A'),('C','G')]
    for x in replacements:
        rev=rev.replace(x[0],x[1])
    for holder in [read,rev]:
        #Search Operation
        lower_main_address,upper_main_address,rel_bias =get_loc(holder)
        i=lower_main_address
        while i<=upper_main_address:
            id = Map[i]-rel_bias # Re-shifting to yield address
            i+=1
            main_address,rel_address=get_Rel_address(id)
            refer = reference[main_address][rel_address:]
            refer=refer[:-1]
            #Extract from blocks
            while len(refer)<len(holder):
                main_address=main_address+1
                refer =refer+ reference[main_address][:-1]
            refer = refer[:len(holder)]

            # Checking the number of mismatch
            flag=True
            if len(holder)!=len(refer):
                flag=False
            else:
                index_break=0
                mistakes=0
                print("Check")
                while index_break<len(holder):
                    if holder[index_break]!=refer[index_break]:
                        mistakes+=1
                    index_break+=1
                    if mistakes>2:
                        flag=False
                        break
            print("Check")
            # Update the list if Flag is High
            if flag==True:
                List_of_matches.append(id)

    # function body - Begins

    positions=List_of_matches

    # function body - Ends
    return positions # list of potential locations at which the read may match the reference sequence.

def WhichExon(positions):
    """
    Input: list of potential locations at which the read may match the reference sequence.
    Output: Update(increment) to the counts of the 12 exons
    IMPORTANT: This function must also handle the following:
        1. cases where the read does not match with the reference sequence
        2. cases where there are more than one matches (match at two exons)
        3. any other special case that you may encounter
    """
    r1,r2,r3,r4,r5,r6,g1,g2,g3,g4,g5,g6 = 0,0,0,0,0,0,0,0,0,0,0,0
    # function body - Begins

    if len(positions)==0: # Edge case
        return np.array([r1,r2,r3,r4,r5,r6,g1,g2,g3,g4,g5,g6])
    
    for id in positions:
        # Checking for Red gene
        if (id>=RedExonPos[0][0] and id<=RedExonPos[0][1]):
            r1=1
        if (id>=RedExonPos[1][0] and id<=RedExonPos[1][1]):
            r2=1
        if (id>=RedExonPos[2][0] and id<=RedExonPos[2][1]):
            r3=1
        if (id>=RedExonPos[3][0] and id<=RedExonPos[3][1]): 
            r4=1
        if (id>=RedExonPos[4][0] and id<=RedExonPos[4][1]): 
            r5=1
        if (id>=RedExonPos[5][0] and id<=RedExonPos[5][1]):
            r6=1


        #Checking Green Gene
        if (id>=GreenExonPos[0][0] and id<=GreenExonPos[0][1]): 
            g1=1
        if (id>=GreenExonPos[1][0] and id<=GreenExonPos[1][1]): 
            g2=1
        if (id>=GreenExonPos[2][0] and id<=GreenExonPos[2][1]):
            g3=1
        if (id>=GreenExonPos[3][0] and id<=GreenExonPos[3][1]): 
            g4=1
        if (id>=GreenExonPos[4][0] and id<=GreenExonPos[4][1]): 
            g5=1
        if (id>=GreenExonPos[5][0] and id<=GreenExonPos[5][1]):
            g6=1
        
        if r1==g1 and g1==1:
            r1,g1=0.5,0.5
        if r2==g2 and g2==1:
            r2,g2=0.5,0.5
        if r3==g3 and g3==1:
            r3,g3=0.5,0.5
        if r4==g4 and g4==1:
            r4,g4=0.5,0.5
        if r5==g5 and g5==1:
            r5,g5=0.5,0.5
        if r6==g6 and g6==1:
            r6,g6=0.5,0.5

      




    # function body - Ends    
    return np.array([r1,r2,r3,r4,r5,r6,g1,g2,g3,g4,g5,g6])


def ComputeProb(ExonMatchCounts):
    """
    Input: The counts for each exon
    Output: Probabilities of each of the four configurations (a list of four real numbers)
    """
    # function body - Begins
    P0,P1,P2,P3=0,0,0,0
    # Method 1
    ## Using generalised KL divergence as metric

    # R_by_G_2=ExonMatchCounts[1]/ExonMatchCounts[2+6]
    # R_by_G_3=ExonMatchCounts[2]/ExonMatchCounts[2+6]
    # R_by_G_4=ExonMatchCounts[3]/ExonMatchCounts[3+6]
    # R_by_G_5=ExonMatchCounts[4]/ExonMatchCounts[4+6]

    # R_G=np.array([R_by_G_2,R_by_G_3,R_by_G_4,R_by_G_5])
    # h1=np.array([0.5,.5,.5,.5])
    # h2=np.array([1,1,0,0])
    # h3=np.array([0.33,.33,1,1])
    # h4=np.array([0.33,.33,.33,1])
    # H=[h1,h2,h3,h4]
    # inv_score=np.zeros(4)
    # for i in range(4):
    #     h_i=H[i]
    #     inv_score[i]=sp.special.kl_div(R_G,h_i).sum()
    # score=np.divide(1,inv_score)
    # P0,P1,P2,P3=sp.special.softmax(score)

    # method 2

    h0=np.array([float(1/3),float(1/3),float(1/3),float(1/3)])
    h1=np.array([float(1/2),float(1/2),0.0,0.0])
    h2=np.array([float(1/4),float(1/4),float(1/2),float(1/2)])
    h3=np.array([float(1/4),float(1/4),float(1/4),float(1/2)])
    Hypothesises={0:h0,1:h1,2:h2,3:h3}

    prob_list = np.zeros(4)
    for index in range(1,5):
        r=ExonMatchCounts[index]
        g=ExonMatchCounts[index+6]
        total=r+g
        for hyp_index,hyp in Hypothesises.items():
            p=hyp[index-1]
            p_power_r=p**r
            q_power_g=(1-p)**g
            prob_list[hyp_index] += float(sp.special.comb(total,r)*p_power_r*q_power_g)
    P0,P1,P2,P3=prob_list/prob_list.sum()



    # function body - ends
    return [P0, P1, P2, P3]




def BestMatch(ListProb):
    """
    Input: Probabilities of each of the four configurations (a list of four real numbers)
    Output: Most likely configuration (an integer). Refer to lecture slides
    """
    # function body - Begins

    MostLikelyConfiguration=None
    MostLikelyConfiguration=np.argmax(ListProb)

    # function body - ends

    return MostLikelyConfiguration # it holds 0, 1, 2, or 3



if __name__ == "__main__":




    # load all the data files
    LastCol = loadLastCol("../data/chrX_last_col.txt") # loads the last column
    RefSeq = loadRefSeq("../data/chrX.fa") # loads the reference sequence
    Reads = loadReads("../data/reads") # loads the reads
    Map = loadMapToRefSeq("../data/chrX_map.txt") # loads the mapping to the reference sequence






    #######   Additional variables #########
    


    #######   Additional variables #########

    bwt= LastCol.split("\n")
    bwt=bwt[:-1]
    bwt=[suit+"\n" for suit in bwt]


    #Generating Milestones (Rank at a position i)

    len_bwt=len(bwt)
    A=np.zeros(shape=len_bwt,dtype=int)
    C=np.zeros(shape=len_bwt,dtype=int)
    G=np.zeros(shape=len_bwt,dtype=int)
    T=np.zeros(shape=len_bwt,dtype=int)

#############Rough####################


    i=1
    A[0]=bwt[0].count('A')
    C[0]=bwt[0].count('C')
    G[0]=bwt[0].count('G')
    T[0]=bwt[0].count('T')

    while i<len_bwt:
        A[i]=A[i-1]+bwt[i].count('A')
        C[i]=C[i-1]+bwt[i].count('C')
        G[i]=G[i-1]+bwt[i].count('G')
        T[i]=T[i-1]+bwt[i].count('T')
        i+=1

    Milestone={'A':A,'C':C,'G':G,'T':T}


    reference=RefSeq.split("\n")
    reference=reference[1:]
    reference=[suit+"\n" for suit in reference]
    reference=reference[:-1]

   
    RefSeq1=RefSeq[6:]
    RefSeq1=RefSeq1.replace("\n","")

    num_of_A,num_of_C,num_of_G,num_of_T=A[-1],C[-1],G[-1],T[-1]

   


    #########################################



    
    # run the functions
    ExonMatchCounts = np.zeros(12) # initialize the counts for exons
    len_of_Reads=len(Reads)
    for i in tqdm(range(len_of_Reads)): # update the counts for exons
        read=Reads[i]
        positions = MatchReadToLoc(read) # get the list of potential match locations
        ExonMatchCounts += WhichExon(positions) # update the counts of exons, if applicable
    ListProb = ComputeProb(ExonMatchCounts) # compute probabilities of each of the four configurations
    MostLikely = BestMatch(ListProb) # find the most likely configuration
    print("Configuration %d is the best match"%MostLikely)
