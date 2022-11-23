class NumMatrix:
    def __init__(self,rows,cols):
        self.mat = []
        for i in range(rows):
            self.mat.append([])
            for _ in range(cols):
                self.mat[i].append(0.0)
    def __getItem__(self,n):
        return self.mat[n]

    def num_rows(self):
        return len(self.mat)
    def num_cols(self):
        return len(self.mat[0])
    def get_val(self,i,j):
        if i > j: return self.mat[i][j]
        else: return self.mat[j][i]
    def set_val(self,i,j,v) :
        if i > j: 
            self.mat[i][j] = v
        else: 
            self.mat[j][i] = v
    def print_mat(self):
        for row in self.mat: print(row)
        print()
    def min_dist_indexes(self):
        m = self.mat[1][0]
        res = (1,0)
        for i in  range(1,self.num_rows()):
            for j in range(i):
                if self.mat[i][j] < m:
                    m = self.mat[i][j]
                    res = (i,j)
        return res

    def add_row(self,newrow): 
        self.mat.append(newrow)

    def add_col(self,newcol): 
        for r in range(self.num_cols()):
            self.mat[r].append(newcol[r])
    def rm_row(self,idx):
        del self.mat[idx]
    def rm_col(self,idx):
        for r in range(self.num_rows()):
            del self.mat[r][idx]
    def copy_mat(self):
        new_mat = NumMatrix(self.num_rows(),self.num_cols())
        for i in range(self.num_rows()):
            for j in range(self.num_cols()):
                new_mat.mat[i][j] = self.mat[i][j]
        return new_mat




def test():
    mat = NumMatrix(3,3)
    print(mat.print_mat())
        
if __name__ == "__main__":
    test()
    
