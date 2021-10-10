
read cat ( LIBPATH, "INTERVALs.mpl" ):
read cat ( LIBPATH, "IntervalPow.mpl" ):
read cat ( LIBPATH, "rootof.mpl" ):
read cat ( LIBPATH, "subres.mpl" ):


SLV := module()
description "Module for operations with real algebraic numbers";

export
    sign_at_1,
    sign_at_2,
    is_leq,
    is_equal,
    is_geq,
    compare,
    is_in_list,
    solve_1,
    solve_2,
    solve_h_2,

    compute_rational_in_between,

    print_list_algnum_1,
    print_list_algnum_2,

    fdisplay_1,
    fdisplay_2;

local
    sign_at_1_internal,
    sgn_intvl,
    is_inside,
    is_rational,
    get_rational,
    compare_rational_rational,
    compare_rootof_rational,
    compare_by_refinement_only,
    compare_by_refinement_only_2,
    copy_rootof;


# Comparison of (a, b)
#    a < b   -1
#    a > b    1
#    a = b    0
    compare_rational_rational := (a::rational, b::rational) -> `if`(a < b, -1, `if`(a > b, 1, 0));

    compare_rootof_rational := proc(alpha, q::rational)
    local sgn_fq;

        if alpha:-is_rational() then
            return compare_rational_rational(alpha:-get_left(), q);
        fi;

        # Maybe q is outside the isolating interval
        # check if q is in the left of J
        if ( q <= alpha:-J[1] )  then return  1; fi;
        # check if q is in the left
        if ( q >= alpha:-J[2] ) then return -1; fi;

        # No q is inside J
        # RECHECK!!!!
        #  f(J[1])  | f(q)  | result
        #     -     |   -   |  +1
        #     -     |   +   |  -1
        #     +     |   -   |  -1
        #     +     |   +   |  +1
        sgn_fq := signum( subs(alpha:-XX = q, alpha:-poly));

        if (sgn_fq = 0) then
            alpha:-display();
            error("We should have computed the rational root already");
        fi;

        if alpha:-get_sgn_left() < 0 then
            if sgn_fq < 0 then return +1;
            else return -1; fi
        else # alpha:-get_sgn_left() > 0
            if sgn_fq < 0 then return -1;
            else return +1; fi
        fi;
    end proc;


    get_rational := proc(a)
        if type(a, rational) then return a; fi;
        return a:-get_midpoint();
    end proc;

    is_rational := proc(a)
        if type(a, rational) or a:-is_rational() then return true; fi;
        return false;
    end proc;

    # TODO: Rewrite to take into account that alpha and beta are
    # are rationals in root_of representation.
    compare := proc(alpha, beta)
    local a, b, Fal, Fbl, Far, Fbr, g, sL, sR;
        if is_rational(alpha) then
            if is_rational(beta) then
                return compare_rational_rational(get_rational(alpha), get_rational(beta));
            fi;
            return -compare_rootof_rational(beta, get_rational(alpha));
        fi;
        if is_rational(beta) then
            return compare_rootof_rational(alpha, get_rational(beta)); fi;

        # Both are algebraic numbers
        # Refine 2 times.
        alpha:-refine(2); beta:-refine(2);


        # If the intervals are disjoint we can decide
        if (alpha:-J[2] <= beta:-J[1]) then return -1 fi;
        if (beta:-J[2] <= alpha:-J[1]) then return 1; fi;

        # The common interval will be [a, b]
        # decide the left endpoint of the common interval;
        if ( alpha:-J[1] >= beta:-J[1] ) then
            a := alpha:-J[1];
            # Maybe beta is in [ Jb[1], a]
            Fbl := signum(eval(beta:-poly(), beta:-XX = a));
            if ( Fbl * beta:-get_sgn_left() <= 0) then return 1; fi;
        else
            a := beta:-J[1];
            Fal := signum(eval(alpha:-poly, alpha:-XX = a));
            if ( Fal * alpha:-get_sgn_left() <= 0) then return -1; fi;
        fi;

        # decide the right endpoint of the isolating interval
        if ( alpha:-J[2] <= beta:-J[2] ) then
            b := alpha:-J[2];
            # print("b", alpha:-J, beta:-J, b);
            Fbr := signum( eval( beta:-poly, beta:-XX = b));
            # print("b", alpha:-J, beta:-J, b, Fbr);
            if ( Fbr * beta:-get_sgn_right() <= 0) then return -1; fi;
        else
            b := beta:-J[2];
            Far := signum(eval(alpha:-poly, alpha:-XX = b));
            if ( Far * alpha:-get_sgn_right() <= 0) then return 1; fi;
        fi;

        if (a >= b) then
            error "something is very wrong!!!!";
        fi;

        # update the isolating intervals
        alpha:-J := [a, b];
        beta:-J := [a, b];

        # So both alpha and beta lie in [a, b]
        #
        if (alpha:-XX <> beta:-XX) then error("the polynomials should have the same variable!") fi;
        g := gcd( alpha:-poly, beta:-poly);
        sL := signum( eval(g, alpha:-XX = a));
        sR := signum( eval(g, alpha:-XX = b));

        if ( sL * sR  < 0 ) then return 0; fi;

        # They are not equal. Refine forever.
        while true do
            alpha:-refine(2);
            beta:-refine(2);

            if ( alpha:-J[2] <= beta:-J[1] ) then return -1; fi;
            if ( alpha:-J[1] >= beta:-J[2] ) then return  1; fi;
        od;
    end proc;

    is_leq   := (alpha, beta) -> is(compare(alpha, beta) <= 0);
    is_equal := (alpha, beta) -> is(compare(alpha, beta)  = 0);
    is_geq   := (alpha, beta) -> is(compare(alpha, beta) >= 0);


    # We assume that the numbers are not equal
    # We use it to sort the distinct real roots of a univariate_poly
    # Output :  true   if alpha <= beta
    #           false  if alpha >= beta
    compare_by_refinement_only := proc(alpha, beta)
        while true do
            if ( alpha:-J[2] <= beta:-J[1] ) then return true; fi;
            if ( alpha:-J[1] >= beta:-J[2] ) then return false; fi;
            alpha:-refine(2);
            beta:-refine(2);
        od;
    end proc;


# notice that a and b are pairs of algebraic numbers
# i.e a = [ alpha_1, alpha_2 ] and b = [ beta_1, beta_2 ]
# We use it to sort the distinct real roots of a bivariate polynomial system
    compare_by_refinement_only_2 := proc(a, b)
        # a[1]:-display();
        while true do
            if ( a[1]:-J[2] < b[1]:-J[1] ) then  return true;  fi;
            if ( a[1]:-J[1] > b[1]:-J[2] ) then  return false; fi;

            if ( a[2]:-J[2] > b[2]:-J[1] ) then  return true;  fi;
            if ( a[2]:-J[1] < b[2]:-J[2] ) then  return false; fi;

            # Refine intervals
            a[1]:-refine(); a[2]:-refine();
            b[1]:-refine(); b[2]:-refine();
        od;
    end proc:

    ## Return true if the algebraic number alpha in the list (of
    ## algebraic numbers) L
    is_in_list := proc(alpha, L)
    local i;
        [seq( SLV:-compare(alpha, L[i]), i=1..nops(L))];
        ListTools:-Search(0, %);
        if (% > 0) then return true; fi;
        return false;
    end proc:


# Input:
#   f : polynomial
#   alpha : rootof
# OUTPUT:
#   sign( f(alpha) )
    sign_at_1 := proc(f_in, alpha)
    local f, c, XX, r;

        if (f_in = 0) then return 0; fi;
        if (degree(f_in) = 0) then return signum(lcoeff(f_in)); fi;

        XX := indets(f_in)[1];
        # print("XX", XX, f_in);
        c := content(f_in, XX);
        f := f_in / c;

        r :=  sign_at_1_internal(f, alpha);
        return r;
    end proc;


    sign_at_1_internal := proc(f, alpha)
    local XX, g, sL, sR, sgn, rr, tpoly, Ia;

        XX := indets(f)[1];
        if alpha:-is_rational() then
            return signum(subs(XX=alpha:-J[1], f));
        fi;

        tpoly := subs(alpha:-XX = XX, alpha:-poly);

        g := gcd(f, tpoly);
        sL := subs(XX = alpha:-J[1], g);
        sR := subs(XX = alpha:-J[2], g);

        if ( sL * sR < 0 ) then
            # In this case alpha is a root of f
            return 0;
        fi;
        sgn := -2;
        to 100 do # Try a FEW times with interval arithmetic.
#            alpha:-refine();
            Ia := INTERVAL( alpha:-J[1] .. alpha:-J[2] );
            rr := evalr(subs(XX=Ia, f));
            #rr := ( evalr( unapply(f, XX)(INTERVAL(
            #alpha:-J[1]..alpha:-J[2]))));
#            t := [op(op(rr))];
#           print("RR", signum(t[1]), signum(t[2]));
            # print("rr", rr);
            sgn := sgn_intvl(rr);
            if sgn = -2 then
                alpha:-refine(100);
#                next;
#                break;
            else
#                printf("finished.\n");
#                printf ("Returning: %d.\n\n", sgn);
                return sgn;
            fi;
        od;
        print("an error here");
    end proc;



# Bivariate sign_at
#
# Assumptions:
#   - f \in Z[x, y]
    sign_at_2 := proc(f, alpha, beta, XX := 'x', YY := 'y')
    local Sq, sz, tl, tr, Vl, Vr, rr, result, vars, Ia, Ib, tpoly, RR, src, i, IIa;

        # Check for the number of variables of f
        vars := indets(f);
        if ( nops(vars) = 0 ) then error "Something very strange happened";  fi;

        # print("vars", vars, f);

        if ( nops(vars) = 1 ) then
            if ( vars[1] = XX ) then
                return sign_at_1(f, alpha);
            fi;
            if ( vars[1] = YY ) then
                return sign_at_1(f, beta);
            fi;
            error "We assume that the bivariate polynomials are \in Z[x, y]";
        fi;
        if nops(vars) >= 3 then error "We need a bivariate polynomial \in Z[x, y]"; fi;


        # Check if alpha and/or beta is a rational number
        if ( type(alpha, rational) ) then
            if ( type( beta, rational) ) then
                print("here");
                return signum( subs({XX = alpha, YY = beta}, f));
            fi;
            return sign_at_1(numer(subs(XX = alpha, f)), beta);
        fi;
        if ( type(beta, rational) ) then
            return sign_at_1(numer(subs(YY = beta, f)), alpha);
        fi;

        # Check if they are rationals represented by an algebraic number
        if ( alpha:-is_rational() ) then
            if ( beta:-is_rational() ) then
                return signum(subs({XX = alpha:-J[1], YY = beta:-J[1]}, f));
            fi;
            return sign_at_1(numer( subs(XX = alpha:-J[1], f)), beta);
        fi;
        if ( beta:-is_rational() ) then
            return sign_at_1(numer( subs(YY = beta:-J[1], f)), alpha);
        fi;



# alpha:-refine(0):
        # beta:-refine(0):
        # alpha:-display(); prinf( " "); beta:-display(); printf( "\n");

        # First check with interval arithmetic and approximations
        to 2 do
            Ia := INTERVAL( alpha:-J[1] .. alpha:-J[2] );
            Ib := INTERVAL(  beta:-J[1] .. beta:-J[2]  );
            #print("Ia, Ib", Ia, Ib);
            rr := evalr(subs(YY=Ib, collect(subs(XX=Ia, f), YY, evalr)));
            # rr := evalr(unapply(f, XX, YY)(Ia, Ib));

            # t := [op(op(rr))];
            # print("RR2", signum(t[1]), signum(t[2]));
            # print("rr2", rr);
            result := sgn_intvl(rr);
            if ( result <> -2 ) then return result; fi;
            alpha:-refine(100);
            beta:-refine(100);
        od;

        #  error("DDDD");
        # alpha:-display(); prinf( " "); beta:-display(); printf( "\n");

        # So, f is a bivariate polynomial
        # and alpha, beta Ran :-)

        # we want tpoly and f to have the same variables.
        tpoly := subs(indets(alpha:-poly)[1] = XX, alpha:-poly);

        RR := RegularChains:-PolynomialRing([XX, YY]);
        src := RegularChains:-ChainTools:-SubresultantChain(f, tpoly,  XX, RR);

        # print(f, tpoly, XX);
        # print("src", ListTools:-Reverse(src[subresultant_chain_vector]));

        Sq := ListTools:-Reverse(src[subresultant_chain_vector]);
        Sq[1] := f;
        Sq[1] := tpoly;
#        Sq := [f, tpoly, seq(RegularChains:-ChainTools:-SubresultantOfIndex(i, src, RR),
#                             i = degree(tpoly, XX)-1..0, -1) ] :
#        print("Sq");
#        Sq := PRS:-StHa(tpoly, f, XX);

        sz := nops(Sq);
        # print("Seq ok");
        # tl := map(numer,  map(apply, map(unapply, Sq, XX),
        # alpha:-J[1]));
        # use the initial interval of alpha, IIa, because it has
        # number of smaller bitsize
        tl := map(numer,  map(apply, map(unapply, Sq, XX), alpha:-sJ[1]));
        Vl := PRS:-Var( map(sign_at_1, tl, beta));

        #print("first leg ok");
        tr := map(numer,  map(apply, map (unapply, Sq, XX), alpha:-sJ[2]));
        Vr := PRS:-Var( map(sign_at_1, tr, beta));

        #print("second leg ok");
        return signum( (Vl - Vr) *  alpha:-get_sgn_deriv() );
    end proc;






# The main function
# This is the only way to construct real algebraic numbers.
#
# CAUTION:
#    The assumption is that f \in Z[x]
    solve_1 := proc(F)
    local sol, m, f, L, g, i, ta, J, lst_polys, lst_mult, lst_intvl, tmp;

        if ( nops( indets( f)) > 1 ) then  error "Not univariate polynomial"; fi;
        if ( nops( indets( f)) = 0 ) then  return []; fi;

        if (lcoeff(F) < 0) then f := numer(-F); else f := numer(F); fi;

        sol := []:

        L:= factors(f)[2];
        lst_polys := map( numer, map2( op, 1, L));
        lst_mult := map2(op, 2, L);

        m := 1:   # the multiplicity;
        for g in lst_polys do
            if degree(g) = 1 then
                ta := solve(g);
                J := [ta, ta]:
                sol := [ op(sol), rootof(g, J, m) ];
            else
                tmp  := RootFinding[Isolate](g, output=interval, digits=1):
                lst_intvl := map2(op, 2, tmp);
                sol := [ op(sol), seq( rootof(g, lst_intvl[i], lst_mult[m]), i=1..nops(lst_intvl)) ];
            fi;

            m := m + 1;
        od:

        return sort(sol, compare_by_refinement_only);
        #        return sort( sol, compareByRefining);
        # return ListTools[Reverse]( sort( sol, compareByRefining));
    end proc;



    solve_2 := proc(F_in, XX := 'x', YY := 'y')
    local sol, sx, sy, F, Lab, Rx, Ry, Fx, Fy, i, j, k, PP, rtx, rr, ix, iy, rty ;

        F := F_in;
        sol := RootFinding[Isolate](F, [XX, YY], digits=10, output=interval);

        # Isolating intervals and boxes
        sol := map( (a) -> [ op(op(a)[1])[2], op(op(a)[2])[2] ], sol);
        sx := map( (a) -> a[1], sol);
        sy := map( (a) -> a[2], sol);
#         print("sol", sol);

        Rx := resultant(F[1], F[2], YY);
        Ry := resultant(F[1], F[2], XX);

        PP := [ seq([0, 0], i=1..nops(sol)) ];

        Fx := map( (f) -> numer(op(f)[1]), factors(Rx)[2]):
        rtx := []:
        for j from 1 to nops(Fx) do
            # Evaluate at the endpoints and obtain the sign
            map( (a)-> [ signum(subs(XX=a[1], Fx[j])), signum(subs(XX=a[2], Fx[j])) ], sx);
            # Check if there are opposite signs at the endpoints
            map( (a) -> a[1]*a[2], %);
            # there is a root if there are opposite signs or the evaluation is 0
            rr := ListTools:-Flatten([ListTools[SearchAll](-1, %),  ListTools[SearchAll](0, %)]);
            for k in rr do PP[k][1] := j; od;
            rtx := [op(rtx), rr ];
        od:

        Fy := map( (f) -> numer(op(f)[1]), factors(Ry)[2]):
        rty := []:
        for j from 1 to nops(Fy) do
            # Evaluate at the endpoints and obtain the sign
            map( (a)-> [ signum(subs(YY=a[1], Fy[j])), signum(subs(YY=a[2], Fy[j])) ], sy);
            # Check if there are opposite signs at the endpoints
            map( (a) -> a[1]*a[2], %);
            # there is a root if there are opposite signs or the evaluation is 0
            rr := ListTools:-Flatten([ListTools[SearchAll](-1, %),  ListTools[SearchAll](0, %)]);
            for k in rr do PP[k][2] := j; od;
            rty := [op(rty), rr];
        od:

        Lab := []:
        for k from 1 to nops(PP) do
            ix := PP[k][1];
            iy := PP[k][2];
            Lab := [ op(Lab),
                     [rootof(Fx[ix], sx[k], 1), rootof(Fy[iy], sy[k], 1)] ];
        od:
        return sort(Lab, compare_by_refinement_only_2);;


#    return [ seq( [ La[PP[i][1]], Lb[PP[i][2]] ], i=1..nops(PP)) ];

    end proc;


    copy_rootof := (a)-> rootof(a:-poly, a:-J, a:-mult);



    # INPUT:
    #  JJ : INTERVAL
    # Output:
    #  the sign of JJ. -2 if JJ contains 0
    sgn_intvl := proc(JJ)
    local J;
        J := [op(op(JJ))];
        # if (nops(J) = 1) then return signum(J); fi;
        if (J[1] = J[2]) then return signum(J[1]); fi;

        if J[1] >= 0 then return 1; fi;
        if J[2] <= 0 then return -1 fi;

        return -2;
    end proc:


    # is interval J1 in interval J2
    is_inside := (J1, J2) -> is(J2[1] <= J1[1]) and is(J1[1] <= J2[2]) and
                             is(J2[1] <= J1[2]) and is(J1[2] <= J2[2]):


    # Real solving of associated polynomials.
    # We need it for PTOPO
    # The system is symmetric, so we output one resultant and a list
    # of pairs of indices, that correspond to solutions of the system.
    solve_h_2 := proc(h, SS := 's', TT := 't')
    local TSsols, Rt, T2sol, Tsol, p, i, tt, L, lpos, j, Ssol, L3, r;

        TSsols := SLV:-solve_2(h, TT, SS):
        Rt := resultant(h[1], h[2], SS):
        Tsol := SLV:-solve_1(Rt):

                # print("start"); for r in Tsol do r:-display(); od; print("end");

        p := [ ]:
        #for i from 2 to 2 do
        for i from 1 to nops(Tsol) do
            tt  := Tsol[i];

            T2sol := [seq( TSsols[i][1], i=1..nops(TSsols) )];
            L := map( SLV:-compare, T2sol, tt);
            lpos := [ ListTools[SearchAll](0, L) ] ;
            # print("lpos", %, tt:-get_approximation(), T2sol[2]:-get_approximation(), SLV:-compare(tt, T2sol[2]));
            if nops(lpos) > 0 then
                for j in lpos do
                    Ssol := TSsols[j][2]; Ssol:-set_var(t);
                    L3 := map( SLV:-compare, Tsol, Ssol);
                    member( 0, L3, 'jy');
                    p := [ op(p), [i, jy]];
                od;
            fi;
        od:

        return p, Tsol;
    end proc:

    # nb : how many rationals in between
    compute_rational_in_between := proc(alpha, beta, nb::integer := 1)
    local r, h;
        # alpha:-refine(30);
        # beta:-refine(30);

        # the number should NOT be equal
        if is_equal(alpha, beta) then
            WARNING("This should never happen. Compute a rational between two equal numbers!");
            return get_rational(alpha);
        fi;

        if type(alpha, rational) then
          if type(beta,rational) then #both rationals
            r := compare(alpha, beta);

            if r <= 0 then
                # alpha <= beta
                h := (beta - alpha) / (2^(nb)+1);
                return seq(alpha + h*k, k=1..2^(nb));
            else
              # beta <= alpha
              h := (alpha - beta) / (2^(nb)+1);
              return seq(beta + h*k, k=1..2^(nb));
           fi;
          else # a rational, b algebraic
          r := compare(alpha, beta);

            if r <= 0  then
                # alpha <= beta
                h := (beta:-J[1] - alpha) / (2^(nb)+1);
                return seq(alpha + h*k, k=1..2^(nb));
                # return alpha:-J[2] + h;
                # return (alpha:-J[2] + beta:-J[1])/2;
              else
                # beta <= alpha
                h := (alpha - beta:-J[2]) / (2^(nb)+1);
                return seq(beta:-J[2] + h*k, k=1..2^(nb));
                #return beta:-J[2] + h;
                #return (alpha:-J[1] + beta:-J[2])/2;
              fi;

          fi:
        else

          if type(beta,rational) then #a algebraic , b rational

            r := compare(alpha, beta);

            if r <=0 then
                # alpha <= beta
                h := (beta - alpha:-J[2]) / (2^(nb)+1);
                return seq(alpha:-J[2] + h*k, k=1..2^(nb));
          else
            # beta <= alpha
            h := (alpha:-J[1] - beta) / (2^(nb)+1);
            return seq(beta + h*k, k=1..2^(nb));
         fi;

          else #both are algebraic numbers

            r := compare_by_refinement_only(alpha, beta);

            if r = true then
                # alpha <= beta
                h := (beta:-J[1] - alpha:-J[2]) / 2^(nb);
                return seq(alpha:-J[2] + h*k, k=1..2^(nb));
                # return alpha:-J[2] + h;
                # return (alpha:-J[2] + beta:-J[1])/2;
            else
              # beta <= alpha
              h := (alpha:-J[1] - beta:-J[2]) / 2^(nb);
              return seq(beta:-J[2] + h*k, k=1..2^(nb));
              #return beta:-J[2] + h;
              #return (alpha:-J[1] + beta:-J[2])/2;
           fi;

          fi:

      fi:
    end proc;

    print_list_algnum_1 := proc(L::list)
    local alpha;
        for alpha in L do alpha:-display(); od:
        return;
    end proc;

    print_list_algnum_2 := proc(L::list)
    local p;
        for p in L do
            p[1]:-display();
            p[2]:-display();
            printf("\n");
        od;
        return;
    end proc;

    fdisplay_1 := proc(L)
    local s;
        for s in L do
            printf("[ %+5.10f ] \n",s:-get_approximation());
        od;
    end proc;

    fdisplay_2 := proc(L)
    local s;
        for s in L do
            printf("[ %+5.10f, %+5.10f ] \n", s[1]:-get_approximation(), s[2]:-get_approximation());
        od;
    end proc;


end module: # SLV
