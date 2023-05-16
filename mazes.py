#! /usr/bin/env python3
''' Run cool maze generating algorithms. '''
import random

class Cell():
    ''' Represents a single cell of a maze.  Cells know their neighbors
        and know if they are linked (connected) to each.  Cells have
        four potential neighbors, in NSEW directions. 
    
    '''
   
    def __init__(self, row, column):
        assert row >= 0
        assert column >= 0
        self.row = row
        self.column = column
        self.links = {}
        self.north = None
        self.south = None
        self.east  = None
        self.west  = None


        
    def link(self, cell, bidirectional=True):
        ''' Carve a connection to another cell (i.e. the maze connects them)'''
        assert isinstance(cell, Cell)
        self.links[cell] = True
        if bidirectional:
            cell.link(self, bidirectional=False)
        
    def unlink(self, cell, bidirectional=True):
        ''' Remove a connection to another cell (i.e. the maze 
            does not connect the two cells)
            
            Argument bidirectional is here so that I can call unlink on either
            of the two cells and both will be unlinked.
        '''
        assert isinstance(cell, Cell)
        del self.links[cell]
        if bidirectional:
            cell.unlink(self, bidirectional=False)
            
    def is_linked(self, cell):
        ''' Test if this cell is connected to another cell.
            
            Returns: True or False
        '''
        assert isinstance(cell, Cell)
        
        #if find cell in  self.links return True, otherwise False
        if cell in self.links:
            return True
        else:
            return False
 
    def all_links(self):
        ''' Return a list of all cells that we are connected to.'''
        return self.links

    
    def link_count(self):
        ''' Return the number of cells that we are connected to.'''
        return len(self.links)
        
    def get_neighbors(self):
        ''' Return a list of all geographical neighboring cells, regardless
            of any connections.  Only returns actual cells, never a None.
              
        '''         
 
        self.neighbors=[]
        if self.north != None:
            self.neighbors.append(self.north)
        if self.south != None:
            self.neighbors.append(self.south)
        if self.east != None:
            self.neighbors.append(self.east)
        if self.west != None:
            self.neighbors.append(self.west)      
      
        return self.neighbors

    
                
    def __str__(self):
        return f'Cell at {self.row}, {self.column}'
        

class Grid:
    ''' A container to hold all the cells in a maze. The grid is a 
        rectangular collection, with equal numbers of columns in each
        row and vis versa.
    '''
    
    def __init__(self, num_rows, num_columns):
        
        assert num_rows > 0
        assert num_columns > 0
        self.num_rows = num_rows
        self.num_columns = num_columns
        self.grid = self.create_cells()
        self.connect_cells()
        
    def create_cells(self):
        ''' Call the cells into being.  Keep track of them in a list
            for each row and a list of all rows (i.e. a 2d list-of-lists).
            
            Do not connect the cells, as their neighbors may not yet have
            been created.
            
        for columns in range(0,self.num_columns):
            for rows in range(0,self.num_rows):
                cell=Cell(rows,columns)
                self.grid[rows][columns]=Cell(rows+1,columns+1)

            
        '''
#       self.grid=[]
        self.grid = [[Cell(i,j) for j in range(self.num_columns)] for i in range(self.num_rows)]              

        return self.grid
                 
                

            
    def connect_cells(self):
        
        ''' Now that all the cells have been created, connect them to 
            each other - set the north/south/east/west attributes.

        '''    

        for column in range(1,self.num_columns-1):
            cell=self.cell_at(0,column)            
            cell.south=self.grid[1][column]
            cell.east=self.grid[0][column+1]
            cell.west=self.grid[0][column-1]

            cell=self.cell_at(self.num_rows-1,column)
            cell.north=self.grid[self.num_rows-2][column]
            cell.east=self.grid[self.num_rows-1][column+1]
            cell.west=self.grid[self.num_rows-1][column-1]

        for row in range(1,self.num_rows-1):
            cell=self.cell_at(row,0)
            cell.south=self.grid[row+1][0]
            cell.east=self.grid[row][1]
            cell.north=self.grid[row-1][0]

            cell=self.cell_at(row,self.num_columns-1)
            cell.south=self.grid[row+1][self.num_columns-1]
            cell.west=self.grid[row][self.num_columns-2]
            cell.north=self.grid[row-1][self.num_columns-1]
        

        for column in range(1,self.num_columns-1):
            for row in range(1,self.num_rows-1):
                cell=self.cell_at(row,column)  
                cell.east=self.grid[row][column+1]
                cell.west=self.grid[row][column-1]
                cell.north=self.grid[row-1][column]
                cell.south=self.grid[row+1][column]
                
        cell=self.grid[0][0]             # set the  connection for 4 cells at the 4 corners
        cell.east=self.grid[0][1]
        cell.south=self.grid[1][0]

        cell=self.grid[self.num_rows-1][0]
        cell.east=self.grid[self.num_rows-1][1]
        cell.north=self.grid[self.num_rows-2][0]                    
        
        cell=self.grid[0][self.num_columns-1]
        cell.south=self.grid[1][self.num_columns-1]
        cell.west=self.grid[0][self.num_columns-2]

        cell=self.grid[self.num_rows-1][self.num_columns-1]
        cell.west=self.grid[self.num_rows-1][self.num_columns-2]
        cell.north=self.grid[self.num_rows-2][self.num_columns-1]
        
                    
        return self.grid
      
        
    def cell_at(self, row, column):
        ''' Retrieve the cell at a particular row/column index.'''
        
        return self.grid[row][column]
        
        
    def deadends(self):
        ''' Return a list of all cells that are deadends (i.e. only link to
            one other cell).
        '''
        deadends_list=[]
        for column in range(0,self.num_columns):
            for row in range(0,self.num_rows):
                cell=self.grid[row][column]
                if len(cell.links) == 1:
                    deadends_list.append(cell)
                    
        return deadends_list
                

                            
    def each_cell(self):
        ''' A generator.  Each time it is called, it will return one of 
            the cells in the grid.
        '''
        for row in range(self.num_rows):
            for col in range(self.num_columns):
                c = self.cell_at(row, col)
                yield c
                
    def each_row(self):
        ''' A row is a list of cells.'''
        for row in self.grid:
            yield row
               
    def random_cell(self):
        ''' Chose one of the cells in an independent, uniform distribution. '''
        randcol = random.randrange(self.num_columns)
        randrow = random.randrange(self.num_rows)
        cell = self.cell_at(randrow, randcol)
        return cell
    
    def size(self):
        ''' How many cells are in the grid? '''
        return self.num_columns * self.num_rows
        
    def set_markup(self, markup):
        ''' Warning: this is a hack.
            Keep track of a markup, for use in representing the grid
            as a string.  It is used in the __str__ function and probably
            shouldn't be used elsewhere.
        '''
        self.markup = markup
        
    def __str__(self):
        ret_val = '+' + '---+' * self.num_columns + '\n'
        for row in self.grid:
            ret_val += '|'
            for cell in row:
                cell_value = self.markup[cell]
                ret_val += '{:^3s}'.format(str(cell_value))
                if not cell.east:
                    ret_val += '|'
                elif cell.east.is_linked(cell):
                    ret_val += ' '
                else:
                    ret_val += '|'
            ret_val += '\n+'
            for cell in row:
                if not cell.south:
                    ret_val += '---+'
                elif cell.south.is_linked(cell):
                    ret_val += '   +'
                else:
                    ret_val += '---+'
            ret_val += '\n'
        return ret_val
        
class Markup:
    ''' A Markup is a way to add data to a grid.  It is associated with
        a particular grid.
        
        In this case, each cell can have a single object associated with it.
        
        Subclasses could have other stuff, of course
    '''
    
    def __init__(self, grid, default=' '):
        self.grid = grid
        self.marks = {}  # Key: cell, Value = some object
        self.default = default
        
    def reset(self):
        self.marks = {}
        
    def __setitem__(self, cell, value):
        self.marks[cell] = value
        
    def __getitem__(self, cell):
        return self.marks.get(cell, self.default)
        
    def set_item_at(self, row, column, value):
        assert row >= 0 and row < self.grid.num_rows
        assert column >= 0 and column < self.grid.num_columns
        cell = self.grid.cell_at(row, column)
        if cell:
            self.marks[cell]=value
        else:
            raise IndexError
    
    def get_item_at(self, row, column):        
        '''Get a cell( at row, column )  martks and return it

        '''
        
        assert row >= 0 and row < self.grid.num_rows
        assert column >= 0 and column < self.grid.num_columns
        cell = self.grid.cell_at(row, column)
        if cell:
            return self.marks.get(cell)
        else:
            raise IndexError
            
    def max(self):
        ''' Return the cell with the largest markup value. '''
        
        return max(self.marks.keys(), key=self.__getitem__)

    def min(self):
        
        ''' Return the cell with the largest markup value. '''
        
        return min(self.marks.keys(), key=self.__getitem__)
    

class DijkstraMarkup(Markup):
    ''' A markup class that will run Djikstra's algorithm and keep track
        of the distance values for each cell.
    '''

    def __init__(self, grid, root_cell, default=0):
        
        ''' Execute the algorithm and store each cell's value in self.marks[]
        '''
        super().__init__(grid, default)
        
        self.grid=grid;
        self.root_cell=root_cell

       
        c=self.root_cell     #  set  cells markup from  root cell
        value=0
        
        cellVisited=[]
        cellVisited.append(c) 
        


        frontier=[]
        frontier.append(c)            # intial frontier list with start_cell
                        
        while len(cellVisited)-1< self.grid.size():
            
            neighbors_unvisited=[]         #  find all cells  in  the frontier
            
            for cell in frontier:
                #cellVisited.append(cell)  #all cells in frontier  are marked as visited
                
                self.set_item_at(cell.row, cell.column, value)  # set the marks for cells in frontier
             
                #find cell's   unvisited linked  neighbors,put them all into neighbors_unvisited
                for neighbor_cell in cell.links:
                    if neighbor_cell not in cellVisited:
                        neighbors_unvisited.append(neighbor_cell)  #??
                cellVisited.append(cell)  #all cells in frontier  are marked as visited
                        
            frontier=[];                         # preparing to move next
            for cell in neighbors_unvisited:    # new frontier with new  unvisited and linked cells
                if cell not in frontier:
                    frontier.append(cell)
            
                
            value +=1    #will move next step and value added by 1
#            print(len(cellVisited), 'how many cell visisted')



            
    def farthest_cell(self):
        ''' Find the cell with the largest markup value, which will
            be the one farthest away from the root_call-self.
            
            Returns: Tuple of (cell, distance)
            
            Return the cell with the largest markup value.
            
        '''
        cell= max(self.marks.keys(), key=self.__getitem__)
        return (cell,self.marks[cell])
    
    

    def shortest_cell(self):

        #return min(self.marks.keys(), key=self.__getitem__)
        
        cell= min(self.marks.keys(), key=self.__getitem__)
        return(cell,self.marks[cell])
    

class ShortestPathMarkup(DijkstraMarkup):
    ''' Given a starting cell and a goal cell, create a Markup that will
        have the shortest path between those two cells marked.  
    '''
    def __init__(self, grid, start_cell, goal_cell, path_marker='*', non_path_marker=' ', end_marker='^'):

            
        super().__init__(grid,start_cell)
        
       
        #self.grid = grid 
        self.start_cell = start_cell
        self.goal_cell = goal_cell
        self.path_marker=path_marker
        self.end_marker=end_marker
        
        
        dm = DijkstraMarkup(self.grid,self.start_cell)

        cellnow=self.goal_cell     

        path=[]     # initial  the path from start cell to goal cell
        path.append(cellnow)     #record  the steps starting from goal_cell

        running=True                
        while running:
            
           #find all linked neighbors
            neighbors_markup={}            # initial a dic for cellnow's all  linked neighbors
            
            for cell in cellnow.links:      #cellnow: the cell  where the foot  is on
                
                mark=dm.marks.get(cell)
                if mark  != None:
                    neighbors_markup.update({cell:mark})   # work out a  dic for  its all linked neighbors
                
                #print( 'neighbor cell and mark', cell, neighbors_markup[cell])
                
            #find the cell in all linked neighbors which  has min markup value
                                    
            next_cell=min(neighbors_markup,key=neighbors_markup.get)
            
            #print('next cell', next_cell )    
               
            path.append(next_cell)    #update the path
                    
            cellnow=next_cell
            #print('cell now',cellnow,'marks',dm.marks[cellnow],'path lenghth',len(path)-1)
            
            if dm.marks.get(cellnow)== 0:   # path reaches the start_cell
                running=False
                
        # markup the path with relative value "*"
        
        for cell in path:
            self.marks[cell]=path_marker
            
        self.marks[self.start_cell]=end_marker
        self.marks[self.goal_cell]=end_marker
   
        
class LongestPathMarkup(ShortestPathMarkup):
    ''' Create a markup with the longest path in the graph marked.
        Note: Shortest path is dependent upon the start and target cells chosen.
              This markup is the longest path to be found _anywhere_ in the maze.
    '''

    def __init__(self, grid, path_marker='*', non_path_marker=' ', end_marker='^'):
       
        start_cell = grid.random_cell()
        
        #print('start cell  ',start_cell)
        
        dm = DijkstraMarkup(grid, start_cell)
        farthest, distance = dm.farthest_cell()
        
        #print('farthest cell from start_cell',  farthest, start_cell, 'distance', distance) 
        
        dm = DijkstraMarkup(grid, farthest)
        next_farthest, distance = dm.farthest_cell()

        #print('next farthest cell from farthest cell',  next_farthest, farthest, 'distance', distance)

        # mark the path from farthest cell to next_farthest cell with relative  markers
        # find  start-end path and mark up the path(shortest path finding)
        
        super().__init__(grid, farthest, next_farthest, path_marker, non_path_marker, end_marker)



class ColorizedMarkup(Markup):
    ''' Markup a maze with various colors.  Each value in the markup is
        an RGB triplet.
    '''

    def __init__(self, grid, channel='R'):
        assert channel in 'RGB'
        super().__init__(grid)
        self.channel = channel
        
    def colorize_dijkstra(self, start_row = None, start_column = None):
        ''' Provide colors for the maze based on their distance from
            some cell.  By default, from the center cell.
        '''
        if not start_row:
            start_row = self.grid.num_rows // 2
        if not start_column:
            start_column = self.grid.num_columns // 2
        start_cell = self.grid.cell_at(start_row, start_column)
        dm = DijkstraMarkup(self.grid, start_cell)
        self.intensity_colorize(dm)
                
    def intensity_colorize(self, markup):
        ''' Given a markup of numeric values, colorize based on
            the relationship to the max numeric value.
        '''
        max = markup.max()
        max_value = markup[max]
        for c in self.grid.each_cell():
            cell_value = markup[c]
            intensity = (max_value - cell_value) / max_value
            dark   = round(255 * intensity)
            bright = round(127 * intensity) + 128
            if self.channel == 'R':
                self.marks[c] = [bright, dark, dark]
            elif self.channel == 'G':
                self.marks[c] = [dark, bright, dark]
            else:
                self.marks[c] = [dark, dark, bright]   
                                       
def binary_tree(self):
    ''' The Binary Tree Algorithm.
      
        This algorithm works by visiting each cell and randomly choosing
        to link it to the cell to the east or the cell to the north.
        If there is no cell to the east, then always link to the north
        If there is no cell to the north, then always link to the east.
        Except if there are no cells to the north or east (in which case
        don't link it to anything.)
      
        
   '''
        
    for column in range(0,self.num_columns):
        for row in range(0,self.num_rows):
            cell=self.cell_at(row,column)

            if cell.north ==None and cell.east !=None:
                cell.link(cell.east)
                
            if cell.east ==None and cell.north !=None:
                cell.link(cell.north)

            if cell.east !=None and cell.north !=None:
                r = random.choice([cell.north,cell.east])
                cell.link(r)                

            
def sidewinder(self, odds=.5):
    ''' The Sidewinder algorithm.
    
        Considers each row, one at a time.
        For each row, start with the cell on the west end and an empty list 
        (the run).  Append the cell to the run list.
        Choose a random number between 0 and 1.  If it is greater 
        than the odds parameter, then add the eastern cell to the run list and
        link it to the current cell.  That eastern cell then becomes the 
        current cell.
        
        If the random number was less than the odds parameter, then you are
        done with the run.  Choose one of the cells in the run and link it to 
        the cell to the north.
        
        Be careful, these instructions don't cover the cases where the row
        is the northernmost one (which will need to be a single, linked run) 
        or for cells at the far east (which automatically close the run)
    '''
    assert odds >= 0.0
    assert odds < 1.0
    
    for columns in range(0,self.num_columns-1):    # most top row linked
        cell=self.grid[0][columns]
        cell.link(cell.east)

        

    for rows in range(1,self.num_rows):
        column = 0                       # initialize  a run while starting a new row
        run = []
        cell = self.cell_at(rows,column)   #cell 
        run.append(cell)
       
        while column < self.num_columns:
            randnum = random.random( )
            if randnum > odds:
                if cell.east != None:
                    run.append(cell.east)
                    column += 1
                    cell.link(cell.east)    # carve to the east neighbored cell
                    cell = cell.east        # walk to the east  neighbored cell                 
                else:
                    cell_picked=random.choice(run)     #end  up  the run since it is end of the row
                    cell_picked_north=self.cell_at(cell_picked.row-1,cell_picked.column)
                    cell_picked.link(cell_picked_north)
                    break
                              
            else:
                                                             #end  up the run
                cell_picked=random.choice(run)                 #random cell to be selected
                cell_picked_north=self.cell_at(cell_picked.row-1,cell_picked.column) # north  cell
                       
                cell_picked.link(cell_picked_north)    # radom cell in the run links with  its northern one
                
                if  cell.east != None:  #initialize the next run
                    column += 1
                    run=[]                             # empity the run list for next run at a row
                    cell=cell.east
                    run.append(cell)
        

          
                
def aldous_broder(self):
    ''' The Aldous-Broder algorithm is a random-walk algorithm.
    
        Start in a random cell.  Choose a random direction.  If the cell
        in that direction has not been visited yet, link the two cells.
        Otherwise, don't link.
        Move to that randomly chosen cell, regardless of whether it was
        linked or not.
        Continue until all cells have been visited.
        
    '''    

    randcell=self.random_cell()
      
    cellVisited = [] 
    cellVisited.append(randcell)
    cellnow = randcell
    iteration_count=0
       
    while len(cellVisited) < self.size():

        iteration_count +=1
        arounds=cellnow.get_neighbors()
        r = random.choice(arounds)
       
        if r not in cellVisited:
            cellnow.link(r)
            cellnow = r
            cellVisited.append(cellnow)
        else:
            cellnow=r
          
    print(f'Aldous-Broder executed on a grid of size {self.size()} in {iteration_count} steps.')
    
    
def wilson(self):
    ''' Wilson's algorithm is a random-walk algorithm.
    
        1) Choose a random cell.  Mark it visited.
        2) Choose a random unvisited cell (note, this will necessarily not be the 
          same cell from step 1).  Perform a "loop-erased" random walk until
          running into a visited cell.  The cells chosen during this random
          walk are not yet marked as visited.
        3) Add the path from step 2 to the maze.  Mark all of the cells as visited.
          Connect all the cells from the path, one to each other, and to the 
          already-visited cell it ran into.
        4) Repeat steps 2 and 3 until all cells are visited.
        
        Great.  But, what is a "loop-erased" random walk?  At each step, one 
        random neighbor gets added to the path (which is kept track
        of in order).  Then, check if the neighbor is already in the path.  If 
        so, then the entire loop is removed from the path.  So, if the 
        path consisted of cells at locations (0,0), (0,1), (0,2), (1,2), (1,3),
        (2,3), (2,2), and the random neighbor is (1,2), then there is a loop.
        Chop the path back to (0,0), (0,1), (0,2), (1,2) and continue 
        
        BTW, it  may be easier to manage a  list of unvisited cells, which 
        makes it simpler to choose a random unvisited cell, for instance.   
    '''
    dest_cell=self.random_cell()    #destination cell
    temp_path=[]              #temp path
    loops_removed=0           # counter for loops removed
    random_choices=0        ## counter for calling random_cell  
     
    cellVisited = []
    cellVisited.append(dest_cell) #dest_cell has been visited
  
    while len(cellVisited) < self.size():   # still has cell to be visited
        
        start_flag=True
        while start_flag:
            cell_now=self.random_cell()      # starting  point (an unvisited cell)
            random_choices +=1
            if cell_now not in cellVisited:
                start_flag=False             # find an unvisited cell as starting cell
        
        temp_path.append(cell_now)
        
         
        path_flag=True
        while path_flag:
            r = random.choice(cell_now.get_neighbors())     #finding current neighbours r
            random_choices +=1
            if r not in cellVisited:
                if r not in temp_path:
                    cell_now = r                         # move to its neighbour
                    temp_path.append(cell_now)           # path becomes longer
                    
                else:                               # erasing the loop path
                    r_index=temp_path.index(r)           # find where r was in the path
                    temp_path=temp_path[:r_index+1]      # ignore cells after r in the path
                    cell_now=r                           # step back to the  cell at the loop crossed  point
                    loops_removed +=1                     #  add 1 to  loops erased 
               
            if r in cellVisited:   # walk is over and all cells in path will be connected and marked Visited
                
                steps_num=len(temp_path)
                
                for n in range(0,steps_num-1):
                    temp_path[n].link(temp_path[n+1])   #link all cells in the path
                    cellVisited.append(temp_path[n])
                temp_path[steps_num-1].link(r)           # link the last cell to the already-visited cell r
                cellVisited.append(temp_path[steps_num-1])  # make the last cell into cellVisited
                temp_path=[]                              # prepare for next path
                path_flag=False                          # this path is over,  now next path                        
                               
     
               
    print(f'Wilson executed on a grid of size {self.size()} with {random_choices}', end='')
    print(f' random cells choosen and {loops_removed} loops removed')
    

def hybird_abplusw(self):
    '''hybrid version of AB + Wilson's that is as fast as you can make it. It should s
        tart with Aldous-Broder, but after some number of runs (how many?) it should
        switch to Wilson's algorithm. Make sure to do some measurements to ensure you've
        made it run faster.
    '''
    randcell=self.random_cell()
      
    cellVisited = [] 
    cellVisited.append(randcell)
    cellnow = randcell
    iteration_count=0

    hybrid_num=int(self.size()/2)             #find the left number which balance AB+Winson  benefits
       
    while len(cellVisited) < self.size()-hybrid_num:

        iteration_count +=1
        arounds=cellnow.get_neighbors()
        r = random.choice(arounds)
       
        if r not in cellVisited:
            cellnow.link(r)
            cellnow = r
            cellVisited.append(cellnow)
        else:
            cellnow=r
          
    print(f'Aldous-Broder  partially({self.size()-hybrid_num} cells by AB) executed on a grid of size {self.size()} in {iteration_count} steps.')

    # above is the AB to find part of maze

    #below is the Wison to find the remaining of maze

    dest_cell=self.random_cell()    #destination cell
    temp_path=[]              #temp path
    loops_removed=0           # counter for loops removed
    random_choices=0        # counter for calling random_cell
    
    
     
    
    cellVisited.append(dest_cell) #dest_cell has been visited
  
    while len(cellVisited) < self.size():   # still has cell to be visited
        
        start_flag=True
        while start_flag:
            cell_now=self.random_cell()      # starting  point (an unvisited cell)
            random_choices +=1
            if cell_now not in cellVisited:
                start_flag=False             # find an unvisited cell as starting cell
        
        temp_path.append(cell_now)
        
         
        path_flag=True
        while path_flag:
            r = random.choice(cell_now.get_neighbors())     #finding current neighbours r
            random_choices +=1
            if r not in cellVisited:
                if r not in temp_path:
                    cell_now = r                         # move to its neighbour
                    temp_path.append(cell_now)           # path becomes longer
                    
                else:                               # erasing the loop path
                    r_index=temp_path.index(r)           # find where r was in the path
                    temp_path=temp_path[:r_index+1]      # ignore cells after r in the path
                    cell_now=r                           # step back to the  cell at the loop crossed  point
                    loops_removed +=1                     #  add 1 to  loops erased 
               
            if r in cellVisited:   # walk is over and all cells in path will be connected and marked Visited
                
                steps_num=len(temp_path)
                
                for n in range(0,steps_num-1):
                    temp_path[n].link(temp_path[n+1])   #link all cells in the path
                    cellVisited.append(temp_path[n])
                temp_path[steps_num-1].link(r)           # link the last cell to the already-visited cell r
                cellVisited.append(temp_path[steps_num-1])  # make the last cell into cellVisited
                temp_path=[]                              # prepare for next path
                path_flag=False                          # this path is over,  now next path                        
                               
     
               
    print(f'Wilson  partially ({hybrid_num} cells by Wison) executed on a grid of size {self.size()} with {random_choices}', end='')
    print(f' random cells choosen and {loops_removed} loops removed')

    

def recursive_backtracker(self, start_cell=None):
        
    '''
    Recursive Backtracker is a high-river maze algorithm.
    
        1) if start_cell is None, choose a random cell for the start
        2) Examine all neighbors and make a list of those that have not been visited
           Note: you can tell it hasn't been visited if it is not linked to any cell
        3) Randomly choose one of the cells from this list.  Link to and move to that 
           neighbor
        3a) If there are no neighbors in the list, then you must backtrack to the last
            cell you visited and repeat.
            
        Suggestion: Use an explicit stack.  You can write this implicitly (in fact,
        the code will be quite short), but for large mazes you will be making lots of 
        function calls and you risk running out of stack space.
    '''
    stack=[]                         # recording the path
    cell_now=self.random_cell()      # starting  point (an unvisited cell)
    stack.append(cell_now)
    
    cellVisited = []
    cellVisited.append(cell_now)

    while  len(cellVisited) !=self.size():
        
        neighbors_unvisited=[]                       #  initial the list and get a list unvisited  of neighbors
        for cell in cell_now.get_neighbors():       # work out a list collecting unvisited neighbours 
            if cell not in cellVisited:
               neighbors_unvisited.append(cell)
        if len(neighbors_unvisited) !=0:                #  there is at least 1 unvisited neighbor
            r = random.choice(neighbors_unvisited)      #   randomly choice one unvisisted  neighbor:r
            stack.append(r)                            #  path become longer
            cell_now.link(r)
            cellVisited.append(r)                       # mark r cell  with visited
            cell_now=r                                  #   step on the new cell
        else:
            stack.pop()                        # remove the last cell from the path
            cell_now=stack[-1]                 # step back to the new last cell
            
  
