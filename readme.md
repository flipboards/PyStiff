## PyStiff

A python interface of FEM simulation. It handles input file (commands) and data file in text format.

## Usage

Test the example file:

    python pystiff.py example/input.txt

### Input Commands

- __read_data__ [filename]: Read problem file;
- __read_material__ [filename]: Read material data file;
- __apply__ [load/bc]: Apply load or boundary condition from problem file;
- __assemble__ [-p problem_type] [-b bc_type]: Assemble the structure.
  - Problem type: Can only be elastic;
  - Boundary condition type: Can be reduce/toones/largenum;
- __solve__: Solve the problem;
- __eval__ [-p x,y]: Evaluate the displacement and stress of one point;

### Problem File

A block file. Block start with _"BEGIN [ITEM]"_, end with _"END"_.

Blocks must exist:

  - __COORD__: Coordination of nodes. Format is
      
        node_name x(,y)(,z)

    node name can be any different string.
  - __ELEM__: Elements. Format is

        node_name_a (node_name_b node_name_c ...)
  
Optional blocks:
  - __BC__: Boundary condition. Format is

        node_name1 (node_name2 ...) x,(y, ...)


      If a coord is replaced by _"x"_, it is ignored. For example, "(x,5)" means a BC equaling 5 on y.
      If a node has been assigned two BCs on the same coordination, only the last one is kept.

  - __LOAD__: Load. Format is
    
        node_name1 (node_name2 ...) x,(y, ...)

     Loads are added linearly on nodes.

### Material File

Format is "name=value". First column should be the correct name shown below, and value must be a number.

Avaliable materials type:
  - Elastic:
      - poisson: Poisson ratio;
      - elastic_modulus
      - thickness
      - density
