#Perform Three-matrix Affine gap global alignment in Python (adapted from Smith-Waterman alignment)
#give initial values of match, mismatch, opening gap penalty and ordinary gap in sequence
seqmatch = 1
seqmismatch= -1
seqgap = -1
gapopen = -5 
#test sequence    
seq_x = "CCTGCCC"
seq_y = "CCAAGCCCC"
class score_matrix:
    def __init__(self, M, Ix, Iy):
        self.M=M
        self.Ix=Ix
        self.Iy=Iy

#create 3 matrix,give the initial values for the first line and colomn
def create_matrix(rows, cols):
    score_matrix.M = [[0 for col in range(cols+1)] for row in range(rows+1)]
    score_matrix.Ix= [[0 for col in range(cols+1)] for row in range(rows+1)]
    score_matrix.Iy= [[0 for col in range(cols+1)] for row in range(rows+1)]
    score_matrix.Ix[0][0] = gapopen
    score_matrix.Iy[0][0] = gapopen
    score_matrix.M[0][0] = 0
    for j in range(1,cols + 1):
        score_matrix.M[0][j] = -float("inf")
        score_matrix.Ix[0][j] = -float("inf") 
        score_matrix.Iy[0][j] = score_matrix.Iy[0][j - 1] + seqgap
    for i in range(1,rows + 1):
        score_matrix.M[i][0] = -float("inf")
        score_matrix.Iy[i][0] = -float("inf") 
        score_matrix.Ix[i][0] = score_matrix.Ix[i -1][0] + seqgap
    
    return score_matrix

#cite from class code
def sc(a,b):
    sc = seqmatch if a == b else seqmismatch
    return sc

#create trace back table
def create_track(rows, cols):
    global track_table
    track_table = [[0 for col in range(cols+1)] for row in range(rows+1)]
    for i in range(rows+1):
        track_table[i][0] = "Ix"        
    for j in range(cols+1):
        track_table[0][j] = "Iy"
    track_table[0][0] = "M"
    return track_table

#calculate the function
def cal_matrix(i,j):
    #if (i,j) is a match or mismatch    
    m_m = score_matrix.M[i-1][j-1]+sc(seq_x[i-1],seq_y[j-1])
    ix_m=score_matrix.Ix[i-1][j-1]+sc(seq_x[i-1],seq_y[j-1])
    iy_m=score_matrix.Iy[i-1][j-1]+sc(seq_x[i-1],seq_y[j-1])
    score_matrix.M[i][j]=max(m_m,ix_m,iy_m)    
    #if (i,j) is an insert in x    
    m_startx = score_matrix.M[i-1][j]+gapopen+seqgap
    ix_extendx=score_matrix.Ix[i-1][j]+seqgap
    score_matrix.Ix[i][j]=max(m_startx,ix_extendx)
    #if (i,j) is an insert in y    
    m_starty = score_matrix.M[i][j-1]+gapopen+seqgap
    iy_extendy=score_matrix.Iy[i][j-1]+seqgap
    score_matrix.Iy[i][j]=max(m_starty,iy_extendy)
    #give trace back value
    decision = max(score_matrix.M[i][j], score_matrix.Ix[i][j], score_matrix.Iy[i][j])
    if (decision == score_matrix.M[i][j]):
        track_table[i][j] = "M" 

    elif (decision == score_matrix.Ix[i][j]):
        track_table[i][j] = "Ix" 
        
    elif (decision == score_matrix.Iy[i][j]):
        track_table[i][j] = "Iy"           

#builds table used for traceback
def build_matrix(cal_matrix,score_matrix):
    rows=len(seq_x)+1
    cols=len(seq_y)+1
    for i in range(1, rows):
        for j in range(1, cols):
            cal_matrix(i,j,score_matrix)
    return score_matrix

#print out the traceback of the best scoring alignment
def print_traceback(seq_x, seq_y):
    align1 = ""
    align2 = ""
    i = len(seq_x) 
    print("\nlength of seq_X:",i)
    j = len(seq_y) 
    print("length of seq_y:", j)
   
    while i > 0 or j > 0:
        if (track_table[i][j] == "M"):
            align1 = seq_x[i - 1] + align1
            align2 = seq_y[j -1] + align2
            i = i - 1
            j = j -1
        elif (track_table[i][j] == "Ix"):
            align1 = seq_x[i - 1] + align1
            align2 = '-' + align2
            i = i -1
        elif (track_table[i][j] == "Iy"):
            align1 = '-' + align1
            align2 = seq_y[j - 1] + align2
            j = j - 1
        else:
            print('error')
            return
    print("alignment results:")
    print(align1)
    print(align2)
    return

#all included to call
def perform_affine_gap():
    score_matrix=create_matrix(len(seq_x), len(seq_y))
    track_table = create_track(len(seq_x), len(seq_y))
    #print the initial table for Match, insertion in seq_x and insertion in seq_y
    print("\nInitial M")
    for i in score_matrix.M:
        print(i)
    print("\nInitial Ix")
    for i in score_matrix.Ix:
        print(i)
    print("\nInitial Iy")
    for i in score_matrix.Iy:
        print(i)

    #calculate the score matrix
    for i in range(1, len(seq_x) + 1):
        for j in range(1, len(seq_y) + 1):
            cal_matrix(i, j)
    #print the calculated score  matrix
    print("\nM")
    for i in score_matrix.M:
        print(i)
    print("\nIx")
    for i in score_matrix.Ix:
        print(i)
    print("\nIy")
    for i in score_matrix.Iy:
        print(i)    
    print("\nTrack")
    for i in track_table:
        print(i)
    print_traceback(seq_x,seq_y)
    print("ended")
    return
perform_affine_gap()
