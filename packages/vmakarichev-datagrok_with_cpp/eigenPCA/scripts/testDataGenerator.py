""" Data generator for testing PCA.
"""

import csv
import random
import math

# Parameters of pseudo-random distributions
UNIFORM_INT_MIN_VALUE = -100
UNIFORM_INT_MAX_VALUE = 100
UNIFORM_FLOAT_MIN_VALUE = -100.00
UNIFORM_FLOAT_MAX_VALUE = 100.00

def generateUniformData(nameOfFile, numOfRows, numOfIntColumns, numOfFloatColumns, precision=2):
    """  Generate pseudo-random columns and save them to file.
         Each column has a uniform distribution.      
    """
    with open(nameOfFile, 'w', newline='') as file:
        writer = csv.writer(file)

        # create head of the table
        headOfTable = []
        for i in range(numOfIntColumns):
            headOfTable.append('ints' + str(i + 1))        
        for f in range(numOfFloatColumns):
            headOfTable.append('floats' + str(f + 1))        
        writer.writerow(headOfTable)

        divisor = 10 ** precision

        # generate main data
        for r in range(numOfRows):
            row = []

            # generate int values
            for i in range(numOfIntColumns):
                row.append(random.randint(UNIFORM_INT_MIN_VALUE, UNIFORM_INT_MAX_VALUE))

            # generate float values
            for f in range(numOfFloatColumns):
                value = random.uniform(UNIFORM_FLOAT_MIN_VALUE, UNIFORM_FLOAT_MAX_VALUE)               
                row.append(float(math.floor(value * divisor)) / divisor)

            writer.writerow(row)
    
def main():
    """
    Simple console interface.          
    """
    try:
        print("Simple generator of colums with uniform distribution.")
        
        name = input("  Enter name of file: ")
        numOfRows = int(input("  Enter num of rows: "))
        numOfIntColumns = int(input("  Enter num of int columns: "))
        numOfFloatColumns = int(input("  Enter num of float columns: "))
        precision = int(input("  Enter precision of float values: "))

        print("Generating...")
        generateUniformData(name, numOfRows, numOfIntColumns, numOfFloatColumns, precision)     
        print("Done!")

    except OSError:
        print("OS ERROR: check paths and file names!")
        

if __name__ == '__main__':
    main()
       