
from fractions import Fraction

#method using absorbing Markov Chains
"""

Input: Matrix. Initial row is initial state. Terminating states are all 0's.
       All indices indicate a weighted probability of sorts that indicate a 
       likelihood to go to a different state. Want to find the probability
       distribution of terminating states.
Output: List of ints, representing this distrubtion. Formatting explained below 
"""
def answer(m):
    #get the terminating states
    numstates = len(m)
    if numstates == 0:
        return []
    if numstates == 1:
        return [1,1]
    sample_term_state = [0 for i in range(numstates)]
    terminating_states = []
    for i in range(numstates):
        if m[i] == sample_term_state:
            terminating_states.append(i)

    ntermstates = len(terminating_states)
    if ntermstates == 0:
        return [1] + [0]*(ntermstates-1) + [1]
    matrix = generate_identity_matrix(ntermstates)#set up the standard matrix
    transition_states = [i for i in range(numstates) if i not in terminating_states]


    for row in matrix:
        row.extend([0]*(numstates-ntermstates)) #add the necessary trailing zeroes
    
    #add the submatrices below
    for i in transition_states:
        row = []
        for index in terminating_states:
            row.append(m[i][index])
        for index in transition_states:
            row.append(m[i][index])
        matrix.append(row)



    #matrix finished, now normalize to fraction, rename it standard

    standard = []
    for row in matrix:
        denom = sum(row)
        temp = [Fraction(elem, denom) if denom != 0 else 0 for elem in row]
        standard.append(temp)



    #build submatrices
    Q = []
    R = []
    for i in range(-1 * len(transition_states), 0):
        Q.append(standard[i][-1*len(transition_states):])
        R.append(standard[i][:-1*len(transition_states)])


    #get the fundamental matrix
    F = getMatrixInverse(matrix_subtraction(generate_identity_matrix(len(transition_states)), Q))

    
    #limiting matrix
    FR = matrix_multiply(F, R)

    print FR
    if len(FR) == 0:
        return []
        
    result = FR[0]

    #format as required: instead of fractions: numerators, followed by denominator as last elem of array
    denoms = [f.denominator for f in result]
    denom = lcm(*denoms)

    result = [f.numerator * (denom/f.denominator) for f in result]
    result.append(sum(result))
    return result



def gcd(*numbers):
    from fractions import gcd
    return reduce(gcd, numbers)


def lcm(*numbers):  
    def lcm(a, b):
        return (a * b) // gcd(a, b)
    return reduce(lcm, numbers, 1)
 
 
def generate_identity_matrix(dim):
    result = []
    for i in range(dim):
        result.append([1 if k == i else 0 for k in range(dim)])
    return result
 
def matrix_subtraction(m1, m2):
    for i in range(len(m1)):
        m1[i] = [m1[i][j] - m2[i][j] for j in range(len(m1[0]))]
    return m1
 

def matrix_multiply(a,b):
    zip_b = zip(*b)
    # uncomment next line if python 3 : 
    # zip_b = list(zip_b)
    return [[sum(ele_a*ele_b for ele_a, ele_b in zip(row_a, col_b)) 
             for col_b in zip_b] for row_a in a]

def transposeMatrix(m):
    t = []
    for r in range(len(m)):
        tRow = []
        for c in range(len(m[r])):
            if c == r:
                tRow.append(m[r][c])
            else:
                tRow.append(m[c][r])
        t.append(tRow)
    return t

def getMatrixMinor(m,i,j):
    return [row[:j] + row[j+1:] for row in (m[:i]+m[i+1:])]

def getMatrixDeterminant(m):
    #base case for 2x2 matrix
    if len(m) == 2:
        return m[0][0]*m[1][1]-m[0][1]*m[1][0]

    determinant = 0
    for c in range(len(m)):
        determinant += ((-1)**c)*m[0][c]*getMatrixDeterminant(getMatrixMinor(m,0,c))
    return determinant

#borrowed from stackPusher (StackOverflow)
def getMatrixInverse(m):
    determinant = getMatrixDeterminant(m)
    #special case for 2x2 matrix:
    if len(m) == 2:
        return [[m[1][1]/determinant, -1*m[0][1]/determinant],
                [-1*m[1][0]/determinant, m[0][0]/determinant]]

    #find matrix of cofactors
    cofactors = []
    for r in range(len(m)):
        cofactorRow = []
        for c in range(len(m)):
            minor = getMatrixMinor(m,r,c)
            cofactorRow.append(((-1)**(r+c)) * getMatrixDeterminant(minor))
        cofactors.append(cofactorRow)
    cofactors = transposeMatrix(cofactors)
    for r in range(len(cofactors)):
        for c in range(len(cofactors)):
            cofactors[r][c] = cofactors[r][c]/determinant
    return cofactors


def main(): #testcase
    print answer([
        [1, 1, 1, 0, 1, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 1, 1, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 1, 1, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 1, 0, 1, 1, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 1, 0, 1, 0, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    ])

if __name__ == '__main__':
    main()