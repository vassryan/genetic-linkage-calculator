#takes data about gene linkage and determines the number of map units (recombination frequency) for three genes
# take the test crosses and the resulting offspring as inputs, print the recombination frequencies and locations of genes (could use some fancy art thing for this)

# use isupper method to check if letter is uppercase
# use split("*") and split()
# sources https://www.ncbi.nlm.nih.gov/books/NBK22073/
# https://passel2.unl.edu/view/lesson/18b30fa2ff29/4

# number at top of file shows number of three point test crosses to analyze
import numpy as np
from numpy.core.numeric import count_nonzero
from numpy.lib import twodim_base
import pdb
# Note: input file would have other data points added if the program were functional
# I kept the amount of data small to make it easier to test the code and debug

input_file = input("Enter the name of the file containing three point test cross data: ")
output_file = input("Enter the name that you want to give the output file: ")


class recombinant:
    def __init__(self, SymbolList,parent,ThreeDArray):
        self.SymbolList = SymbolList
        self.parent = parent
        self.ThreeDArray = ThreeDArray
    def SingleRecomb(self):
        
        # compare to first parent
        for i in range((len(self.SymbolList))):
            for j in range(3):
                if self.SymbolList[i][j] != self.parent[0][j]:
                    self.ThreeDArray[i][j][0] =1 
                else:
                    continue
        # compare to second parent
        for i in range((len(self.SymbolList))):
            for j in range(3):
                if self.SymbolList[i][j] != self.parent[1][j]:
                    self.ThreeDArray[i][j][1] = 1
                else:
                    continue
        
    def FirstSecondDist(self, ValueList, total):
        # Split 3D array into 2 2D arrays, array 0 is parent 1, array 1 is parent 2
        TwoDArray = np.dsplit(self.ThreeDArray,2)
        Parent1_array = TwoDArray[0]
        Parent2_array = TwoDArray[1]

        # make list that stores the indices for the letters that were different
        Indices_of_Difference_Parent1 = []
        Indices_of_Difference_Parent2 = []
        for i in range(6):
            # count number of differences between group and parent1/2
            parent1_counts = np.count_nonzero(Parent1_array[i]==1)
            parent2_counts = np.count_nonzero(Parent2_array[i]==1)
            
            # Maybe could add index of single recombinants to list when 1st two letters are different from first two parent letters
            
            SingleRecombFirstTwo = []
            if parent1_counts ==2:
                # need to do more research on how to properly use this method, the indexing specifically is very odd
                y = np.where(Parent1_array[i] == 1)
                if y[0,:,0] == 1 and y[0,:,1] == 1:
                    SingleRecombFirstTwo.append(i)
            print(SingleRecombFirstTwo)

            # np.where should find index of the letter that was different
            if parent1_counts == 2 and parent2_counts == 1:
                x = np.where(Parent1_array[i] == 1)
                # this part is not giving correct values
                x_list = [x[:][0][0],x[:][0][1]]
                Indices_of_Difference_Parent1.append(x_list)
            elif parent2_counts == 2 and parent1_counts == 1:
                x2 = np.where(Parent2_array[i] == 1)
                x_list2 = [x2[:][0][0],x2[:][0][1]]
                Indices_of_Difference_Parent2.append(x_list2)
            # initialize list that will help me determine the index of the groups that differed from either parent1 or 2 by 2
            # if parent1 differed by 2, then add 1 to the list, otherwise add 2, can later take index of the 1/2 to see during which loop of the outer for loop it was added
            # 
            # val_list_index = []
            # if parent1_counts == 2:
            #     val_list_index.append(1)
            # elif parent2_counts == 2:
            #     val_list_index.append(2)
            # # List containing the indices to use in value list
            # ValueIndices_parent1 =[]
            # ValueIndices_parent2 =[]

            # for h in val_list_index:
            #     if h == 1:
            #         ValueIndices_parent1.append(val_list_index.index(h))
            #     else:
            #         ValueIndices_parent2.append(val_list_index.index(h))
                            
        print(Indices_of_Difference_Parent1)
        print(Indices_of_Difference_Parent2)

        # sum will be the total number of offspring that demonstrated the FirstSecond recombination event
        FirstTwoSum = 0
        FirstThirdSum = 0
        # Iterate through Indicies_of_Difference list to add appropriate values from value_list to the sum
        # for j in range(len(Indices_of_Difference_Parent1)):
        #     if (Indices_of_Difference_Parent1[j][0] == 0 and Indices_of_Difference_Parent1[j][1] == 2) \
        #         or (Indices_of_Difference_Parent1[j][0] == 2 and Indices_of_Difference_Parent1[j][1] == 0):
        #         FirstThirdSum += ValueList[ValueIndices_parent1[j]]
        #     elif (Indices_of_Difference_Parent1[j][0] == 0 and Indices_of_Difference_Parent1[j][1] == 1) \
        #         (Indices_of_Difference_Parent1[j][0] == 1 and Indices_of_Difference_Parent1[j][1] == 0):
        #         FirstTwoSum += ValueList[ValueIndices_parent1[j]]
        # for k in range(len(Indices_of_Difference_Parent2)):
        #     if (Indices_of_Difference_Parent2[k][0] == 0 and Indices_of_Difference_Parent2[k][1] == 2) \
        #         or (Indices_of_Difference_Parent2[k][0] == 2 and Indices_of_Difference_Parent2[k][1] == 0):
        #         FirstThirdSum += ValueList[ValueIndices_parent2[k]]
        #     elif (Indices_of_Difference_Parent2[k][0] == 0 and Indices_of_Difference_Parent2[k][1] == 1) \
        #         (Indices_of_Difference_Parent2[k][0] == 1 and Indices_of_Difference_Parent2[k][1] == 0):
        #         FirstTwoSum += ValueList[ValueIndices_parent2[k]]


        recomb_frequency_first_two = FirstTwoSum/total
        recomb_frequency_first_third = FirstThirdSum/total
        print(recomb_frequency_first_third)
        print(recomb_frequency_first_two)
    


infile = open(input_file, "r")
outfile = open(output_file,"w")

test_cases = int(infile.readline())

for i in range(test_cases):
        # getting the letters that represent the different genes
    symbols = infile.readline()
    symbols = symbols.split()
    first_symbol = symbols[0]
    second_symbol = symbols[1]
    third_symbol= symbols[2]
    symbol_list = []
    value_list =[]
    for line in range(8):
        # need way to determine when next set of data begins
        if line == "\n":
            continue
        # splitting up the symbols and the value where value is index 3
        # symbols go in symbol_list, values go in value_list
        data = infile.readline()
        group = data.split("*")
        symbol_list.append(" ".join(group[:3]))
        # add values to the value_list, removes \n if present, otherwise just adds value
        # converts string into integer for use in calculations
        if "\n" in group[3]:
            value_list.append(int(group[3][:-1]))
        else:
            value_list.append(int(group[3]))

    # Dictionary that includes keys which are symbols and values which are offspring numbers
    group_dictionary = {symbol_list[0]: value_list[0],symbol_list[1]:value_list[1],symbol_list[2]:value_list[2],
        symbol_list[3]:value_list[3],symbol_list[4]:value_list[4],symbol_list[5]:value_list[5],symbol_list[6]:value_list[6],symbol_list[7]:value_list[7]
    
    }
    # total number offspring in the experiment
    total_offspring = value_list[0] + value_list[1] +value_list[2]+value_list[3]+value_list[4]+value_list[5]+value_list[6]+value_list[7]
    # list containing the symbols for the two groups of offpsring wtih parental genotypes
    parent_list = [symbol_list[0].split(), symbol_list[1].split()]
    
    # Getting list of only the single recombinants
    
    symbol_list = [i.split() for i in symbol_list]
    # Getting list of only the double recombinants
    symbol_list_no_parent = symbol_list[2:8]
  
    # Setting up 3D array to determine which genes are different from parents
    ThreeDimArray = np.zeros((6,3,2))
    Recomb = recombinant(symbol_list_no_parent,parent_list,ThreeDimArray)
    # Make sure to not include double recombinants in calculations until end
    Recomb.SingleRecomb()
    Recomb.FirstSecondDist(value_list,total_offspring)
    # IF CODE WOULD HAVE WORKED, WOULD'VE CHANGED FILE PARAMETER IN PRINT STATEMENTS TO BE EQUAL TO OUTPUT FILE

    # recomb_first_two = Recomb.FirstSecondDist(value_list,total_offspring)* 100
    # wild_type_symbol_list = [first_symbol, second_symbol, third_symbol]
    # print("Recombination Frequency for {0} and {1}: {2:0.2f}%".format(wild_type_symbol_list[0],wild_type_symbol_list[1],recomb_first_two))



infile.close()
outfile.close()
# At end in output file, could print out a line and show the relative positions of the genes
# based on the data, could have each - in the line equal 1 map unit
# ex: v--------cv---------------ct could represent 8 map units between v and cv and 15 between cv and ct
# Would also want to print # of map units so the viewer doesn't have to count each dash

