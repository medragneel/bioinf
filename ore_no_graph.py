





class Ore_no_graph():
    def __init__(self,g={}) :
        self.graph = g
    # graph fn
    def print_graph(self):
        for v in self.graph.keys():
            print(v,"=> ",self.graph[v])
    def get_nodes(self):
        return list(self.graph.keys())
    def get_edges(self):
        edges = []
        for v in self.graph.keys():
            for u in self.graph[v]:
                edges.append((v,u))
        return edges
    def size(self):
        return len(self.get_nodes()), len(self.get_edges())
    def add_vertx(self,v):
        if v not in self.graph.keys():
            self.graph[v] = [] 
    def add_edge(self,o,d):
        if o not in self.graph.keys():
            self.add_vertx(o)
        if d not in self.graph.keys():
            self.add_vertx(d)
        if d not in self.graph[o]:
            self.graph[o].append(d)
    def get_successors(self,v):
        return list(self.graph[v])
    def get_predecessors(self,v):
        res = []
        for k in self.graph.keys():
            if v in self.graph[k]:
                res.append(k)
        return res
    def get_adj(self,v):
        suc = self.get_successors(v)
        pre = self.get_predecessors(v) 
        res = pre
        for p in suc:
            if p not in res: res.append(p)
        return res
    # degree
    def out_degree(self,v):
        return len(self.graph[v])
    def in_degree(self,v):
        return len(self.get_predecessors(v))
        
    def degree(self,v): 
        return len(self.get_adj(v))
    def all_degrees(self,deg_type="inout"):
        # deg_type = "in" or deg_type = "out" or deg_type="inout"
        degs = {}
        for v in self.graph.keys():
            if (deg_type == "out" or deg_type == "inout"):
                degs[v] = len(self.graph[v])
            else:
                degs[v] = 0
        if (deg_type == "in" or deg_type == "inout"):
            for v in self.graph.keys():
                for d in self.graph[v]:
                    if deg_type == "in" or v not in self.graph[d]:
                        degs[d] = degs[d] + 1
        return degs
    def highest_degrees(self,all_deg=None,deg_type="inout",top=10):
        if all_deg is None:
            all_deg = self.all_degrees(deg_type)
        ord_deg = sorted(list(all_deg.items()),key=lambda x:x[1],reverse=True)
        return list(map(lambda x:x[0], ord_deg[:top]))
    # bfs and dfs search
    def reachable_bfs(self,v):
        l=[v]
        res=[]
        while len(l) > 0:
            node = l.pop(0)
            if node != v: res.append(node)
            for elem in self.graph[node]: 
                if elem not in res and elem not in l and elem != node:
                    l.append(elem)
        return res

    def reachable_dfs(self,v):
        l = [v]
        res=[]
        while len(l) > 0:
            node = l.pop(0)
            if node != v: res.append(node)
            s=0
            for elem in self.graph[node]: 
                if elem not in res and elem not in l :
                    l.insert(s,elem)
                    s += 1

        return res
    def distance(self,s,d):
        if s == d: return 0
        l = [(s,0)] 
        visited = [s]
        while len(l) > 0:
            node,dist = l.pop(0)
            for elem in self.graph[node]:
                if (elem == d ): return dist + 1
                elif elem not in visited:
                    l.append((elem, dist+1))
                    visited.append(elem)
        return None

    def shortest_path(self,s,d):
        if s == d: return 0
        l = [(s,[])] 
        visited = [s]
        while len(l) > 0:
            node,preds = l.pop(0)
            for elem in self.graph[node]:
                if (elem == d ): return preds + [node,elem] 
                elif elem not in visited:
                    l.append((elem, preds +[node]))
                    visited.append(elem)
        return None
    def reachable_with_dist(self,s):
        res=[]
        l = [(s,0)]
        while len(l) > 0:
            node,dist = l.pop(0)
            if node != s:res.append((node,dist))
            for elem in self.graph[node]:
                if not is_in_tuple_list(l,elem) and not is_in_tuple_list(res,elem): 
                    l.append((elem,dist+1))
        return res
    # cycles
    def node_has_cycle(self,v):
        l=[v]
        res = False
        visited = [v]
        while len(l) > 0 :
            node = l.pop(0)
            for elem in self.graph[node]:
                if elem == v: return True
                elif elem not in visited:
                    l.append(elem)
                    visited.append(elem)
        return res
    
    def has_cycle(self):
        res = False
        for v in self.graph.keys():
            if self.node_has_cycle(v): return True
        return res





def is_in_tuple_list(tl, val):
    res = False
    for (x,y) in tl:
        if val == x: return True
    return res

                







def test():
    gr = Ore_no_graph()
    gr.add_vertx(1)
    gr.add_vertx(2)
    gr.add_vertx(3)
    gr.add_vertx(4)
    gr.add_edge(1,2)
    gr.add_edge(2,3)
    gr.add_edge(3,2)
    gr.add_edge(3,4)
    gr.add_edge(4,2)
    print(gr.graph)
    gr.print_graph()
    print(gr.size())
    print (gr.get_successors(2))
    print (gr.get_predecessors(2))
    print (gr.get_adj(2))
    
    print (gr.in_degree(2))
    print (gr.out_degree(2))
    print (gr.degree(2))
    
    print(gr.all_degrees("inout"))
    print(gr.all_degrees("in"))
    print(gr.all_degrees("out"))

    gr2 = Ore_no_graph({1:[2,3,4], 2:[5,6],3:[6,8],4:[8],5:[7],6:[],7:[],8:[]})
    print(gr2.reachable_bfs(1))
    print()
    print(gr2.reachable_dfs(1))
    print()

    print(gr2.distance(1,7))
    print(gr2.shortest_path(1,7))
    print(gr2.distance(1,8))
    print(gr2.shortest_path(1,8))
    print(gr2.distance(6,1))
    print(gr2.shortest_path(6,1))
    
    print(gr2.reachable_with_dist(1))
    print()    
    print(gr.has_cycle())
    print(gr2.has_cycle())


    
    

if __name__ == "__main__":
    test()
