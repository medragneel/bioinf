class Binarytree():
    def __init__(self,val,d=0,left=None,right=None):
        self.val = val
        self.d = d
        self.left = left
        self.right = right

    def print_tree(self):
        self.print_tree_rec(0,"Root")
    def print_tree_rec(self,level,side):
        tabs=""
        for _ in range(level):tabs += "\t"
        if self.val >= 0:
            print(tabs,side,"_value:",self.val)
        else:
            print(tabs,side,"_dist:",self.d)
            if(self.left != None):
                self.left.print_tree_rec(level +1 ,"Left")
                
            if(self.right != None):
                self.right.print_tree_rec(level +1 ,"Right")
    def get_cluster(self):
        res=[]
        if self.val >=0:
            res.append(self.val)
        else:
            if(self.left != None):
                res.extend(self.left.get_cluster())
                
            if(self.right != None):
                res.extend(self.right.get_cluster())
        return res



def test():
    a = Binarytree(1)
    b = Binarytree(2)
    c = Binarytree(3)
    d = Binarytree(4)
    e = Binarytree(-1,2.0,b,c)
    f = Binarytree(-1,1.5,d,a)
    g = Binarytree(-1,4.5,e,f)
    g.print_tree()
    print(f.get_cluster())
    print(g.get_cluster())


if __name__ == "__main__":
    test()
    
                
    
