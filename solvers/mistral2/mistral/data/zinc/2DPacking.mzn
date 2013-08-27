% RUNS ON zinc_fdic_mznlib
% RUNS ON zinc_mip_hl
% RUNS ON mzn_fd
% RUNS ON mzn_mer_fd
% RUNS ON mzn_cd_fd
% RUNS ON mzn_mer_lp
% RUNS ON minizinc_fd
%------------------------------------------------------------------------------%
% 2DPacking.mzn
% Jakob Puchinger
% October 29 2007
% vim: ft=zinc ts=4 sw=4 et tw=0
%------------------------------------------------------------------------------%

% annotations commented out for the moment,
% since minizinc does not know about strings
% 
% annotation colgen_definition(int: S, array[1..S] of int: K,
%           array[1..S] of string: SPNames, array[1..S] of string: SPSolvers);
%
% annotation colgen_orig_var(string SPName);
%
% annotation colgen_subproblem_constraint(string: sp_id);
%
% annotation lp_branch_and_bound(array[int] of var int: vars, string: lp_solver)
%     string: variable_choice, string: split);
%
% annotation colgen_dw(string: master_solver);
%

%------------------------------------------------------------------------------%
% Two-dimensional Bin Packing 
%
% n rectangular items with given height and width have to be packed into
% rectangular bins of size W x H
% It is assumed that the items are sorted according to non-increasing height.


    % Upper Bound on the number of Bins used
int: K;
    % Width of the Bin
int: W;
    % Height of the Bin
int: H;
    % Number of items to be packed
int: N;
    % Widths of the items
array[1..N] of int: ItemWidth;
    % Heights of the items
array[1..N] of int: ItemHeight;

    % variable indicating if bin k is used
array[1..K] of var 0..1: bin :: is_output; % :: colgen_orig_var("sp1"); 

    % variable indicating if item i is in bin k
array[1..K, 1..N] of var 0..1: item :: is_output; % :: colgen_orig_var("sp1"); 

    % The total number of bins has to be minimized
solve minimize 
    sum( k in 1..K ) ( bin[k] ) ;
    % :: colgen_definition(1, [1: K], [1: "sp1"], [1: "mip"])
    % :: lp_branch_and_bound( [ x | x in { bin[k] | k in 1..K } union
    %                                    { item[k, j] | k in 1..K, j in 1..N} ],
    %                         "colgen_dw(cplex)", "most_frac", "std_split" ); 

    % each item has to be packed
constraint
    forall( j in 1..N ) (
        sum( k in 1..K ) ( item[k, j] ) = 1
    );

    % subproblem constraints
constraint
    forall( k in 1..K ) (
%        is_feasible_packing(bin[k], [ j: item[k, j] | j in 1..N ])
        is_feasible_packing(bin[k], [ item[k, j] | j in 1..N ])
            %:: colgen_subproblem_constraint("mip", "sp1", k)
    ); 

    % This predicate defines a feasible packing
    % as a 2-stage guillotine pattern (level pattern).
    % A Bin consists of one or several levels 
    % and each level consist of one or several items.
    % The height of a level is given by its highest (first) item.
    % Variable x[i, j] indicates if item i is put in level j
    % x[j, j] also indicate that level j is used.
    %
    % Note, k is locally fixed, so this is the pattern for the k'th bin. 
    %
predicate is_feasible_packing(var 0..1: s_bin,
        array[1..N] of var 0..1: s_item) = (
    let{array[1..N, 1..N] of var 0..1: x} in (
        forall ( i in 1..N ) (
                % Width of items on level can not exceed W
            sum( j in i..N ) ( ItemWidth[j] * x[i, j] ) <= W * x[i, i]
            /\
                % first item in level is highest
                % XXX do not need to generate those variables (use default)
            forall( j in 1..i-1 ) ( x[i, j] = 0 ) 
        )
        /\    
            % Height of levels on bin can not exceed H
        sum( i in 1..N ) ( ItemHeight[i] * x[i, i] ) <= s_bin * H 
        /\
            % for all items associate item on bin with item on level.
        forall(j in 1..N) (
            s_item[1] = sum( i in 1..j ) ( x[i, j] )
        )
    )
);

    % required for showing the objective function
var int: obj;
constraint
    obj = sum( k in 1..K )( bin[k] );

output 
[ "Number of used bins = ",  show( obj ), "\n"] ++ 
[ "Items in bins = \n\t"] ++ 
[ show(item[k, j]) ++ if j = N then "\n\t" else " " endif |
    k in 1..K, j in 1..N ] ++
[ "\n" ];

%------------------------------------------------------------------------------%
%
% Test data (Not in separate file, so that mzn2lp can handle it).
%

N = 4;
W = 5;
H = 10;
K = 2;
ItemWidth =  [ 1, 1, 2, 3 ];
ItemHeight = [ 4, 4, 4, 3 ];

%N = 6;
%W = 5;
%H = 10;
%K = N;
%ItemWidth =  [ 4, 4, 1, 1, 2, 3 ];
%ItemHeight = [ 6, 5, 4, 4, 4, 3 ];
