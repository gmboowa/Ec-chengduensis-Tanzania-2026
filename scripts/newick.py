"""Minimal Newick parser -> nested nodes with name, length, support, children."""
class Node:
    __slots__=('name','length','support','children','parent','x','y')
    def __init__(self):
        self.name=''; self.length=0.0; self.support=None
        self.children=[]; self.parent=None; self.x=0.0; self.y=0.0
    @property
    def is_tip(self): return not self.children

def parse(text):
    s=text.strip().rstrip(';'); i=0
    def skip():
        nonlocal i
        while i<len(s) and s[i] in ' \t\r\n': i+=1
    def read_node():
        nonlocal i
        n=Node(); skip()
        if i<len(s) and s[i]=='(':
            i+=1
            while True:
                c=read_node(); c.parent=n; n.children.append(c); skip()
                if i<len(s) and s[i]==',': i+=1; continue
                if i<len(s) and s[i]==')': i+=1; break
                raise ValueError('bad newick at %d'%i)
        j=i
        while i<len(s) and s[i] not in ',():': i+=1
        lab=s[j:i].strip().strip("'\"")
        if n.children:
            n.support=float(lab) if lab else None
        else:
            n.name=lab
        if i<len(s) and s[i]==':':
            i+=1; j=i
            while i<len(s) and (s[i].isdigit() or s[i] in '.eE+-'): i+=1
            n.length=float(s[j:i])
        return n
    root=read_node()
    return root

def tips(n):
    if n.is_tip: yield n
    for c in n.children: yield from tips(c)

def ladderize(n, ascending=True):
    if n.is_tip: return 1
    sizes=[(ladderize(c,ascending),c) for c in n.children]
    sizes.sort(key=lambda t:t[0], reverse=not ascending)
    n.children=[c for _,c in sizes]
    return sum(k for k,_ in sizes)

def layout(root):
    """x = distance from root, y assigned to tips in current order."""
    def setx(n, acc):
        n.x=acc+n.length
        for c in n.children: setx(c, n.x)
    root.length=0.0; setx(root, 0.0)
    counter=[0]
    def sety(n):
        if n.is_tip:
            n.y=counter[0]; counter[0]+=1; return n.y
        ys=[sety(c) for c in n.children]
        n.y=(min(ys)+max(ys))/2.0
        return n.y
    sety(root)
    return counter[0]
