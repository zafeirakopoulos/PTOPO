
PTOPO := module()
description "Tools for studying the topology of parametric curves";
option package;

######### Exported variables and procedures #########
export
### Procedures

    # Initialization
    parametrization,

    # Properness
    is_proper,
    make_proper,

    # Compute points of interest
    is_there_point_at_infinity,
    compute_topology,

    ## Various
    get_phi,
    get_assoc,
    get_box,

    # Main actions
    topology,
    draw,

    # Printing
    isolated_points,
    point_at_infinity,
    poles,
    extreme_points,
    cusps,
    double_points,
    boundary_points,
    point_summary,

    # Examples
    example
    ;



######### Local variables and procedures #########
local
    ### Submodules
    Box2D,
    Box3D,

    ### Procedures


    # Reparametrizations
    reparametrize,
    random_reparametrization,
    initialize,

    # Compute points of interest
    compute_poles,
    compute_point_at_infinity,
    compute_extreme_points,
    compute_intersections_with_box_2,
    compute_double_points_and_cusps,
    compute_points,
    compute_isolated_points,

    # Auxiliary
    get_points,
    flatten_to_pairs,
    get_point,
    get_segment,
    print_list_points,
    compute_bounding_box,


    ### Variables
    T, S,        # names of the variables. Default TT = 't' and SS = 's'
    n,           # The ambient dimension of the curve
    NP, DP,       # List of the numerators and denominators polynomials
    phi,         # The rational map phi
    h,           # The associated polynomials
    p_inf,       # point at infinity
    box,
    Ht,           # t values for the solutions of the h-system

    pL ,           # The list of points defining the topology
    Bdt,          # t values for intersection with the bounding box
    Pt,           # t values for  poles
    Dt,           # t values for double points
    Ct,           # t values for cusps
    Et,           # t values for extreme points
    Sq,           # The list of Segments

    IP_lst,       # Isolated points. List of (x, y) coords
    Tt,

    # point indices
    idxSt,        # Categorize to points
    idxSPt,        # Categorize to points
    idxDt,         # indices of double points
    idxCt       # indices of cusps
    ;

    ##############################
    ###### Box2D Submodule #######
    ##############################

    Box2D := proc(il := -1, ir := 1, ib := -1, it := 1)
    description "2D box";
        module()
        export
            left,
            right,
            top,
            bottom;

        local init;

        option
            load = init;

            init := proc()
                left := il;
                right := ir;
                bottom := ib;
                top := it;
            end proc;

            init();
        end module: ## Box
    end proc:

    ##############################
    ###### Box3D Submodule #######
    ##############################

    Box3D := proc(il := -1, ir := 1, ib := -1, it := 1, ifr := 1, ibk := -1)
    description "3D box";
        module()
        export
            left,
            right,
            top,
            bottom,
            front,
            back;

        local init;

        option
            load = init;

            init := proc()
                left := il;
                right := ir;
                bottom := ib;
                top := it;
                front := ifr;
                back := ibk;
            end proc;

            init();
        end module: ## Box
    end proc:


###############################################################################
######### Procedures ##########################################################
###############################################################################

##############################
####### Initialization #######
##############################

# Set the parametrization given a list of rational functions
    parametrization := overload(
        [
            proc( p_input::list, q_input::list, iTT := 't', iSS :='s' ) option overload;
            local p, q;
                p := p_input;
                q := q_input;
                initialize(p, q, iTT, iSS);
            end proc,
            proc( phi_input::list, iTT := 't', iSS :='s') option overload;
            local p, q;
                p := map(numer, phi_input);
                q := map(denom, phi_input);
                initialize(p, q, iTT, iSS);
            end proc,
            proc( ) option overload;
                print(NP);
                print(DP);
                return phi;
            end proc
        ]
                               ):

    # Initializes T, S, NP, ND, phi, n

    initialize := proc(p1::list, q1::list, iTT := 't', iSS :='s')
    local i, L, p, q;
        # Set the variable names
        T := iTT; S := iSS;

        ## Avoid a singular at infinity
        ## TODO We should check if there is a point at infinity, and
        ## then we should reparametrize
        if is_there_point_at_infinity(p1, q1) then
            printf("\t \t  There is a point at infinity. We reparametrize.... \n");
            printf("\n \n");
            p, q := op(map( normal, [random_reparametrization(p1, q1)]));
            # print("pq", p, q);
        else
            p := p1; q := q1;
        fi;

        # Set the ambient dimension of the curve
        n := nops(p);

        # Set phi
        phi := normal([ seq(p[i]/q[i], i=1..n)]):


        for i from 1 to n do
        if is(phi[i], rational)then
            error("***  One of the coordinates is a constant...***");
        fi;
        od:

        # Set numerator and denominator polynomials
        NP := map(numer, phi);
        DP := map(denom, phi);

        # Compute the associated h polynomials, for i \in [n]:
        # h_i(s, t) = ( p_i(t) q_i(s) - p_i(s) q_i(t) ) / (s -t)
        h := numer( simplify( [ seq(
            (NP[i]*subs(T=S, DP[i]) - subs(T=S, NP[i])*DP[i])/(S-T),
            i=1..n) ])):

        if is_proper() then
            printf("\n ----****---- The parametrization is  proper  ----****---- \n");
        else
            printf("\n ----****---- The parametrization is NOT proper  ----****---- \n");
            L := make_proper();
            printf("        The new parametrization is %a \n", L[1]);
            parametrization(L[1], iTT, iSS);
        fi;
    end proc;

    is_there_point_at_infinity := proc(p::list, q::list, iTT := 't')
        return convert([seq( is(degree(p[i], iTT) <= degree(q[i], iTT)), i=1..nops(p))], `and`);
    end proc;


##############################
##### Reparametrizations #####
##############################

    # Will perform the necessary reparametrization in order to have
    # a proper parametrization without singular points at infinity.
    reparametrize := proc(p_input,q_input)
    local p,q;

        p := p_input;
        q := q_input;
        p,q := make_proper(p,q):
        p,q := random_reparametrization(p,q):
        return p,q;
    end proc;

    random_reparametrization := proc(p_input,q_input)
    local t0,p,q;
        t0 := rand(1..15)();
        p := subs(T = (t0*T+1)/(T-t0), p_input):
        q := subs(T = (t0*T+1)/(T-t0), q_input):
        return p,q;
    end proc;


##############################
####### Compute Points #######
##############################
# compute_points
# compute_poles
# compute_point_at_infinity
# compute_double_points_and_cusps
# compute_extreme_points


    compute_points := proc()
    local tL,i;

        # Compute the points of interest
        compute_poles();
        compute_point_at_infinity();
        compute_double_points_and_cusps();
        compute_extreme_points();
        compute_intersections_with_box_2();

        # Compute and sort the parameters for all the points of interest
        tL := ListTools:-MakeUnique( ListTools:-Flatten( [idxCt, idxDt]));

        # List of parameters of all the points,
        # Poles (Pt), Boundary (Bdt), Extreme (Et), special points
        # (indices in tL)

        Tt := [ op(Pt), op(Bdt), op(Et), seq( Ht[i], i in tL) ];
        Tt := sort(Tt, (a, b) -> `if`( is(SLV:-compare(a, b) = -1), true, false) );

        compute_isolated_points();
    end proc;


    # NUMERICAL procedure for computing the isolated points.
    compute_isolated_points := proc()
    local H, tD, TT, tt, ev, Im_lst, eps;
        tD := Digits:
        Digits:= 30:
        H := PTOPO:-get_assoc();
        TT := ListTools:-MakeUnique([ fsolve( resultant(H[1], H[2], s), complex) ]):

        IP_lst := []:
        eps := 10^(-10):
        for tt in TT do
            if ( not is(tt, real) ) then
                ev := subs(t = tt, phi);
                Im_lst := map(Im, ev);
                if (abs(Im_lst[1]) < eps) and (abs(Im_lst[2]) < eps) then
                    IP_lst := [ op(IP_lst), map(Re, ev)];
                fi;
            fi;
        od:
        IP_lst := ListTools:-MakeUnique(IP_lst);
        Digits := tD:
        return IP_lst;
    end proc;

    # Computes the poles and set Pt to contain the corresponding parameter values.
    compute_poles := proc()
    local i,p;
        printf(" * Computing the poles... ");

        # Set Pt to contain the parameter values corresponding to poles.
        Pt := SLV:-solve_1( expand(product(DP[i], i=1..n)));

        # Set the label of the points to "p" for pole.
        for p in Pt do p:-other := "p"; od;

        printf("                           -- END\n");
    end proc;



    # Computes point at infinity and set p_inf to contain the corresponding parameter value.
    compute_point_at_infinity := proc()
    local i,res;
        printf(" * Computing point at infinity... ");

        res := convert([seq( is(degree(NP[i], T) <= degree(DP[i], T)), i=1..n)], `and`);
        if res then p_inf := map(limit, phi, T=infinity); fi;
        printf("                   -- END\n");
    end proc;

    # Requires to know the poles, i.e., Pt.
    compute_double_points_and_cusps := proc()
    local i,p, is_a_pole, alpha, res, p2, j, PP, tc, td;


        printf(" * Computing cusps and double points... ");
        PP, Ht := SLV:-solve_h_2(h):

        ## TODO: CHECK
        # SLV:-fdisplay_1(Ht);

        # Keep only one of [a, b] and [b, a]
        p2 := ListTools:-MakeUnique(PP, 1, (x, y)-> is(x[1]=y[2]) and is(x[2]=y[1]));

        # We should NOT consider the poles for singular (special)  points.
        idxSt := []:
        for j from 1 to nops(p2) do
            is_a_pole := false;
            alpha := Ht[p2[j][1]];
            for p in Pt do
                res := SLV:-compare(p, alpha);
                if  res = 0 then
                    is_a_pole := true;
                    break;
                fi;
            od;
            if (not is_a_pole) then idxSt := [ op(idxSt), p2[j]]; fi;
        od;

        # Categorize wrt to the first coordinate, thus to points.
        # for example [ [[1,2], [1,3], [2, 3] [1,1] ] is a triple
        # point and cusp
        idxSPt := [ ListTools:-Categorize((x,y) -> is(x[1]=y[1]), idxSt) ];
        ##  print("idxSPt", idxSPt);
        idxSPt := [ ListTools:-Categorize((x,y) ->
                                          is(nops(
                                              {op((ListTools:-Flatten(x)))}
                                              intersect
                                              {op((ListTools:-Flatten(y)))}) >= 1), idxSPt)
                  ];
        idxSPt := map(flatten_to_pairs, idxSPt);

        idxDt := ListTools:-MakeUnique(ListTools:-Flatten( select( (x)-> is(x[1] <> x[2]), idxSt)));
        idxCt := ListTools:-MakeUnique(ListTools:-Flatten( select( (x)-> is(x[1] = x[2]), idxSt)));
        Ct := [ seq(Ht[i], i in idxCt) ];
        Dt := [ seq(Ht[i], i in idxDt) ];

        # each entry in idxSPt is a list of indices that corresponds
        # to parameters, that in turn correspond to a point p. For
        # example p could be the list [ [[1,2], [1,3], [2, 3] [1,1] ]
        # which corresponds to a point that is a triple point and cusp
        for p in idxSPt do
            tc, td := selectremove( (x)-> is(x[1] = x[2]), p);
            # parameters for the cusps
            tc := ListTools:-MakeUnique( ListTools:-Flatten( tc));
            # parameters for the double points
            td := ListTools:-MakeUnique( ListTools:-Flatten( td));

            for i in tc do
                Ht[i]:-other := cat( Ht[i]:-other, "c");
            od;
            for i in td do
                Ht[i]:-other := cat( Ht[i]:-other, "d");
            od;
        od;
        printf("             -- END\n");

    end proc;

    # Set Et to the list of parameters for extreme points
    # Requires Pt, Ct, Dt
    compute_extreme_points := proc()
    local Htt, REt, e;
        printf(" * Computing extreme points... ");

        Htt := factor(subs(S=T, h));
        Et := map(SLV:-solve_1, Htt);

        Et := ListTools:-Flatten(Et);
        for e in Et do e:-other := "e"; od;

        Et := ListTools:-MakeUnique( Et, 1, (a, b) -> SLV:-is_equal(a, b));
        Et :=  sort(Et, (a, b) -> `if`( is(SLV:-compare(a, b) = -1), true, false) );
        # SLV:-fdisplay_1(Dt);

        # exclude the poles and the cusps
        REt := []:
        for e in Et do
            if not SLV:-is_in_list(e, Pt)
            and not SLV:-is_in_list(e, Ct)
            and not SLV:-is_in_list(e, Dt) then
                REt := [op(REt), e ];
            fi;
        od;

        Et := REt;
        printf("                      -- END\n");
    end proc;


##############################
####### Bounding Box #########
##############################
# compute_bounding_box
# compute_intersections_with_box_2

    compute_bounding_box := proc()
    local e,r, i, LSP, LSPx, LSPy, LSPz, mx, Mx,
    Rational_parameters, my, My, mz, Mz, a;

        printf(" * Computing intersections with bounding box 2D... ");

        LSP := []:
        Rational_parameters :=[]:
        #Consider rationals in between poles
        if (nops(Pt) > 0) then
          LSP := [subs( T = Pt[1]:-get_midpoint()-1, phi)]: # at a rational before 1rst pole
          Rational_parameters := [Pt[1]:-get_midpoint()-1]:
          i:=1:
          while (i < nops(Pt)) do
            a := SLV:-compute_rational_in_between(Pt[i], Pt[i+1])[1]:
            LSP := [ op(LSP), simplify(subs( T = a, phi))]:
            Rational_parameters := [op(Rational_parameters), a]:
            i := i+1:
          od:
          LSP := [op(LSP), subs( T= Pt[i]:-get_midpoint()+1, phi)]: # add a rational after last pole
          Rational_parameters := [op(Rational_parameters), Pt[i]:-get_midpoint()+1 ]:
        fi:

        # consider all extreme and singular points
        LSP := [ op(LSP), seq( subs(T = e:-get_midpoint(), phi), e in Et),
           seq( subs(T = Ht[i[1][1]]:-get_midpoint(), phi), i in idxSPt) ];


          LSPx := [seq(LSP[i][1], i=1..nops(LSP)) ];
          mx := min(LSPx);
          Mx := max(LSPx);
          LSPy := [seq(LSP[i][2], i=1..nops(LSP))];
          my := min(LSPy);
          My := max(LSPy);
          r := 1;

        if n = 2 then
            box := Box2D( floor(mx)-r, ceil(Mx)+r, floor(my)-r, ceil(My)+r );
        elif n = 3 then
            LSPz := [seq(LSP[i][2], i=1..nops(LSP))];
            mz := min(LSPz);
            Mz := max(LSPz);
            box := Box3D( floor(mx)-r, ceil(Mx)+r, floor(my)-r, ceil(My)+r, floor(mz)-r, ceil(Mz)+r );
        end if;

        printf("  -- END\n");

    end proc;


    # compute the intersections with the bounding box
    compute_intersections_with_box_2 := proc()
    local pol, sq, CLt, Lt, CRt, Rt, CBt, Bt, CTt, Tt, fleft, fright, fbot, ftop, t;

        compute_bounding_box();

        # print("B", box:-left, box:-right, box:-top, box:-bottom);
        pol := NP[1] - box:-left* DP[1];
        CLt := SLV:-solve_1(pol);
        Lt := [];
        ftop := NP[2] - box:-top * DP[2];
        fbot := NP[2] - box:-bottom * DP[2];

        for t in CLt do
            t:-other := "b";
            sq := SLV:-sign_at_1(DP[2], t);
            if ((SLV:-sign_at_1(ftop, t)*sq < 0) and
                (SLV:-sign_at_1(fbot, t)*sq > 0)) then
                Lt := [ op(Lt), t ];
            fi;
        od;


        pol := NP[1] - box:-right* DP[1];
        CRt := SLV:-solve_1(pol);
        Rt := [];
        ftop := NP[2] - box:-top * DP[2];
        fbot := NP[2] - box:-bottom * DP[2];
        for t in CRt do
            t:-other := "b";
            sq := SLV:-sign_at_1(DP[2], t);
            if ((SLV:-sign_at_1(ftop, t)*sq < 0) and
                (SLV:-sign_at_1(fbot, t)*sq > 0)) then
                Rt := [ op(Rt), t ];
            fi;
        od;

        pol := NP[2] - box:-top* DP[2];
        CTt := SLV:-solve_1(pol);
        Tt := [];
        fleft :=  NP[1] - box:-left * DP[1];
        fright := NP[1] - box:-right * DP[1];
        for t in CTt do
            t:-other := "b";
            sq := SLV:-sign_at_1(DP[1], t);
            if ((SLV:-sign_at_1(fleft, t)*sq > 0) and
                (SLV:-sign_at_1(fright, t)*sq < 0)) then
                Tt := [ op(Tt), t ];
            fi;
        od;

        pol := NP[2] - box:-bottom* DP[2];
        CBt := SLV:-solve_1(pol);
        Bt := [];
        fleft :=  NP[1] - box:-left * DP[1];
        fright := NP[1] - box:-right * DP[1];
        for t in CBt do
            t:-other := "b";
            sq := SLV:-sign_at_1(DP[1], t);
            if ((SLV:-sign_at_1(fleft, t)*sq > 0) and
                (SLV:-sign_at_1(fright, t)*sq < 0)) then
                Bt := [ op(Bt), t ];
            fi;
        od;

        Bdt := ListTools:-Flatten([Lt, Rt, Bt, Tt]);
        Bdt :=  sort(Bdt, (a, b) -> `if`( is(SLV:-compare(a, b) = -1), true, false) );
    end proc;

##############################
######### Printing ###########
##############################
# point_summary
# isolated_points
# point_at_infinity
# poles
# boundary_points
# extreme_points
# cusps
# double_points

    print_list_points := proc(L, str)
    local t, coords;
        for t in L do
            coords := subs(T = t:-get_midpoint(), phi);
            printf("%s : [ %+10.6f, %+10.6f],   t: %+10.6f   %s\n",
                   str,
                   evalf(coords[1]), evalf(coords[2]),
                   t:-get_approximation(), t:-other);
        od;
    end proc;


    poles := proc()
    local t;
        printf("* Poles:  \n");
        for t in Pt do
            printf("%s :  t: %+10.6f   %s\n",
                   "",
                   t:-get_approximation(), t:-other);
        od;
    end proc;

    point_at_infinity := proc()
        if is_there_point_at_infinity(NP, DP) then
            printf("* Point at infinity:\n");
            printf("%s : [ %+10.6f, %+10.6f]\n", "", evalf(p_inf[1]), evalf(p_inf[2]));
        fi;
    end proc;

    extreme_points := proc()
        printf("* Extreme points:\n");
        print_list_points(Et,"");
    end proc;

    boundary_points := proc()
        printf("* Boundary points:\n");
        print_list_points(Bdt,"");
    end proc;

    cusps := proc()
        printf("* Cusps:\n");
        print_list_points(Ct,"");
    end proc;

    double_points := proc()
        printf("* Double points:\n");
        print_list_points(Dt,"");
    end proc;

    isolated_points := proc()
    local p;
        printf("* Isolated points:\n");
        for p in IP_lst do
            printf("%s : [ %+10.6f, %+10.6f]    i\n",
                   "",
                   p[1], p[2]);
        od;
        #print_list_points(Dt,"");
    end proc;

    point_summary := proc()
        isolated_points();
        point_at_infinity();
        poles();
        extreme_points();
        boundary_points();
        cusps();
        double_points();
    end proc;

##############################
########## Topology  #########
##############################
# topology
# draw


    # Compute the topology.
    # Sq contains segments that are topologically boring.
    topology := proc()
    local s,tSq, nb1, nb2;

        printf("* Compute topology... \n");
        compute_points();

        # Split to subsequences using the poles
        tSq := [ ListTools:-Split((a) -> is(a:-other="p"), Tt)];
        Sq := ListTools:-SelectFirst(nops(tSq),  (a)-> is(nops(a) > 0), tSq);

        # Check if we should merge the first and the last subsequence
        if nops(Sq) > 1 then
            [seq(s:-other, s in Sq[1])];
            nb1 := nops(ListTools:-SelectFirst(nops(Sq[1]), (a) -> is(a = "b"), %));
            # print(%);
            [seq(s:-other, s in Sq[-1])];
            nb2 := nops(ListTools:-SelectFirst(nops(Sq[-1]), (a) -> is(a = "b"), %));
            if (nb1 = 1 and nb2 = 1) then
                Sq[1] := [ op(Sq[1]), op(Sq[-1]) ];
                Sq[1] := sort(Sq[1], (a, b) -> `if`( is(SLV:-compare(a, b) = -1), true, false) );
                Sq[-1] := []:
                Sq := ListTools:-SelectFirst(nops(Sq),  (a)-> is(nops(a) > 0), Sq);
            fi;
        fi;
        printf("-- END topology\n");
        return Sq;
    end proc;

    draw := overload(
        [
            proc( p_input::list, q_input::list, iTT := 't', iSS :='s', NN::integer := 3 )
            option overload;
                parametrization(p_input, q_input, iTT, iSS):
                topology():
                return draw(NN):
            end proc,

            proc( phi_input::list, iTT := 't', iSS :='s', NN::integer := 3)
            option overload;
                parametrization(phi_input, iTT, iSS):
                topology():
                return draw(NN):
            end proc,

            proc(NN::integer := 2)
            option overload;
            local i, j, points, p, L, r, lst1, lst2;

                pL := []:
                # Now Sq is a list of lists.  We treat each list, L, in Sq,
                # independently.

                printf("Computing points in branches \n");
                # print("Sq", Sq);
                for L in Sq do
                    for i from 1 to nops(L)-1 do
                        if L[i]:-other = "b" and  L[i+1]:-other = "b" then next; fi;
                        # L[i]:-display();
                        points := get_points(L[i], L[i+1], NN);


                        for j from 1 to nops(points)-1 do
                            pL := [op(pL),
                                   get_segment(points[j], points[j+1]),
                                   get_point(points[j]),
                                   get_point(points[j+1])
                                  ];
                        od;
                    od;

                    # We might need to connect the first and the last point
                    if L[1]:-other <> "b" or L[-1]:-other <> "b" then

                        lst1 := [ seq( L[1]:-J[1] - (i), i = 1..NN) ];
                        pL := [ op(pL),
                                get_segment(L[1], lst1[1]),
                                seq(get_segment(lst1[i], lst1[i+1]), i=1..NN-1),

                                seq(get_point(lst1[i]), i=1..NN)
                              ];

                        lst2 := [ seq( L[-1]:-J[2] + i, i = 1..NN) ];
                        pL := [ op(pL),
                                get_segment(L[-1], lst2[1]),
                                seq( get_segment(lst2[i], lst2[i+1]), i=1..NN-1),
                                #        get_segment(lst[1], L[-1]),

                                seq( get_point(lst2[i]), i=1..NN)
                              ] ;

                        pL := [op(pL), get_segment(lst1[NN], lst2[NN]) ];
                        # r := L[-1]:-J[2] + 1;
                        # pL := [op(pL),
                        #        get_segment(L[1], r), get_segment(r, L[-1]),
                        #        get_point(r)
                        #       ];

                    fi;
                od;
                # Isolated points
                for p in IP_lst do
                    pL := [op(pL),
                           plottools:-point(p, color = "red", symbolsize = 9,
                                            symbol = "asterisk") ];
                od;
                printf("Graph computed\n");
                return plots:-display(pL, view=[box:-left .. box:-right, box:-bottom..box:-top], axes=boxed):
            end proc
        ]
                    ):



##############################
########## Auxiliary #########
##############################

    ## Help function to assign the parameters to points
    flatten_to_pairs := proc(lst)
    local i, T1;
        T1 := ListTools:-Flatten(lst);
        [ seq( [T1[2*i+1], T1[2*i+2]], i=0..nops(%)/2 - 1) ];
    end proc:


    get_point := proc(t)
    local w, coords;
        if type(t, rational) then
            w := "o";
            coords :=  subs( T = t, phi);
        else
            w := t:-other;
            coords :=  subs( T = t:-get_midpoint(), phi);
        fi;
        if is(w = "b") or is(w = "o") then
            return plottools:-point(coords, color = "black", symbolsize = 9,
                                    symbol = "asterisk");
        elif is(w = "e") then
            return plottools:-point(coords, color = "khaki",  symbolsize = 9,
                                    symbol = "solidcircle");
        elif is(w = "r") then
            return plottools:-point(coords, color = "blue", symbolsize = 2,
                                    symbol="solidcircle");
        else
            return plottools:-point(coords, color = "red", symbolsize = 11,
                                    symbol = "solidbox");
        fi;
    end proc;


    get_segment := proc(t1, t2)
    local coords, s1, s2;
        if type(t1, rational) then s1 := t1; else s1 := t1:-get_midpoint(); fi;
        if type(t2, rational) then s2 := t2; else s2 := t2:-get_midpoint(); fi;

        coords := [subs(T = s1, phi), subs(T = s2, phi)];
        return plottools:-curve(coords, color = "green", symbolsize = 5);
    end proc;


    get_points := proc(L, R, d)
    local M, r, x, lst, K;
        if (d = 0) then return [L,R]; fi;
        lst := [SLV:-compute_rational_in_between(L, R, d)];
        return [L, op(lst), R];
#    print("lst", lst);
        K := [L]:
        for r in lst do
            M := SLV:-solve_1( denom(r)*x-numer(r))[1];
            M:-other := "r";
            K := [op(K),  M ];
#        print("K", K);
        od;
        K := [op(K), R];
#    print("K", K, nops(K));
        return [L, op(lst), R];

        #r := SLV:-compute_rational_in_between(L,R,2);
        #M:= SLV:-solve_1( denom(r)*x-numer(r))[1];
        #M:-other := "r";
        #return ListTools:-MakeUnique(ListTools:-Flatten([get_points(L,M,d-1),get_points(M,R,d-1)]));

    end proc;



    is_proper := proc()

    if n = 2 then
     return is(degree(gcd(h[1], h[2])) <= 0);
    elif n = 3 then
     return is(degree(gcd(gcd(h[1], h[2]), h[3]) <= 0));
    end if;

    end proc;

    # # Algorithm by Sonia Perez-Diaz to deduce a proper
    # # reparametrization of a rational curve.
    make_proper := proc()
    local TH, m, CoeffList, R, Q, Ck, Cl,G, r, L1, L2, L3, F1, F2, F3, K1, K2, K3,
        found, idx, i, j, np1, np2, np3, nq1, nq2, nq3;

        if n = 2 then
            TH := collect(gcd(h[1], h[2]) * (S-T), S, factor);
        elif n = 3 then
            TH := collect(gcd(gcd(h[1], h[2]), h[3]) * (S-T), S, factor);
        end if;
        m := degree(TH, S);

        if (m = 1) then
            # printf("\n ----****---- The parametrization is already proper  ----****---- \n");
            return [phi, 1];
        fi;

        # how to choose Ck and Cl ?
        CoeffList := [ seq( coeff(TH, S, i), i=0..m) ];
        found := false;
        idx := []:
        for i from 1 to nops(CoeffList)-1 do
            for j from i+1 to nops(CoeffList) do
                if degree( gcd(CoeffList[i], CoeffList[j])) = 0 then
                    # print(CoeffList[i], CoeffList[j]);
                    idx := [i, j]:
                    found := true;
                    break;
                fi;
            od;
            if found = true then break; fi;
        od;

        Ck := CoeffList[idx[1]]; Cl := CoeffList[idx[2]];
        G := S*Cl - Ck;
        r := max(degree(Ck, T), degree(Cl, T));
        R := Ck/Cl;

        F1 := x * DP[1] - NP[1];
        L1 := collect(primpart(resultant(F1, G, T)), x, factor);

        F2 := x*DP[2] - NP[2];
        L2 :=  collect((primpart(resultant(F2, G, T))), x, factor);

        K1 := -coeff(L1, x, m-1)/coeff(L1, x, m)/m;
        np1 := numer(%);
        nq1 := denom(%%);

        K2 := -coeff(L2, x, m-1)/coeff(L2, x, m)/m;
        np2 := numer(%);
        nq2 := denom(%%);

        if n = 2 then
            Q := subs(S = T, [np1/nq1, np2/nq2]);
        elif n = 3 then
            F3 := x*DP[3] - NP[3];
            L3 :=  collect((primpart(resultant(F3, G, T))), x, factor);
            K3 := -coeff(L3, x, m-1)/coeff(L3, x, m)/m;
            np3 := numer(%);
            nq3 := denom(%%);
            Q := subs(S = T, [np1/nq1, np2/nq2, np3/nq3]);
        end if;

        return [Q, normal(R) ];
    end proc;


    get_phi := ()-> phi;
    get_assoc := ()-> h;
    get_box := ()-> box;

#### Examples
    example := proc(n::integer)
    local curves;
        curves := [
            [ [t^5*(2*t - 1), (2*t-1)^3], [t^7, t^7] ],
            [ [ t^2 + 1, 1] ,[t^4 + 1, t^3]],
            [[ (t^2 -1)*(t^4 -1+9*t^2), -(t^8 -2*t^6 +2*t^2 -1-54*t^4) ], [9*t^2*(t^2 + 1), 27*(t^2 + 1)*t^3 ]],
            [ [44+37*t^3 -23*t^2 +87*t, 95-61*t^3 -8*t^2 -29*t],[10+29*t^3 +98*t^2 -23*t, 40+11*t^3 -49*t^2 -47*t ]],
            [[ (t^2 -1)*(t^4 -1+9*t^2), -(t^8 -2*t^6 +2*t^2 -1-54*t^4) ], [9*t^2*(t^2 + 1), 27*(t^2 + 1)*t^3 ]],
            [ [3*t^4+ 4*t^3 + 32*t^2 + 28*t + 99, (t^2 + t + 7)^3],[(t^2+t+7)*(t^2 + 1), (t+6)*(t^2 +1)^2]]
                  ];
        return op(curves[n]);
    end proc;
end module: # PTopo
