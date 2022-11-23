from binary_tree import Binarytree
from numat import NumMatrix

class H_clustering():
    def __init__(self,mat_dists):
       self.mat_dists = mat_dists
    def execute_clustering(self):
        trees = []
        for i in range(self.mat_dists.num_rows()):
            tree = Binarytree(i)
            trees.append(tree)
        table_dist = self.mat_dists.copy_mat()
        for k in range(self.mat_dists.num_rows(),1,-1):
            mins = table_dist.min_dist_indexes()
            i,j = mins[0],mins[1]
            n = Binarytree(-1,table_dist.get_val(i,j)/2.0,trees[i],trees[j])
            if k >2:
                ti = trees.pop(i)
                tj = trees.pop(j)
                dists = []
                for x in range(table_dist.num_rows()):
                    if x != i and x != j:
                        si= len(ti.get_cluster())
                        sj= len(tj.get_cluster())
                        d = (si * table_dist.get_val(i,x) + sj * table_dist.get_val(j,x) / (si + sj))
                        dists.append(d)
                table_dist.rm_row(i)
                table_dist.rm_row(i)
                table_dist.rm_col(j)
                table_dist.rm_col(j)
                table_dist.add_row(dists)
                table_dist.add_col([0] * ( len(dists) + 1 ))
                trees.append(n)
            else: return n


def test():
    m = NumMatrix(5,5)
    m.set_val(0, 1, 2)
    m.set_val(0, 2, 5)
    m.set_val(0, 3, 7)
    m.set_val(0, 4, 9)
    m.set_val(1, 2, 4)
    m.set_val(1, 3, 6)
    m.set_val(1, 4, 7)
    m.set_val(2, 3, 4)
    m.set_val(2, 4, 6)
    m.set_val(3, 4, 3)
    hc = H_clustering(m)
    arv = hc.execute_clustering()
    arv.print_tree()

if __name__ == "__main__":
    test()
