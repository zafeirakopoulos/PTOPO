
# The module of algebraic numbers
rootof := proc(polynomial, interval, multiplicity)
description "A data structure for real algebraic numbers";

    module()
    export
        XX,                # The variable of poly.
        poly,              # Univariate polynomial.
        J, 	               # Exact isolating interval.
        sJ,                # The first isolating interval assigned. It
                           # is never refined. For hacks and speedups.
        mult,              # Multiplicity in original polynomial.

        get_poly,
        get_multiplicity,
        get_degree,
        get_approximation,
        get_interval,
        get_left, get_right,
        get_midpoint,
        get_sgn_left,
        get_sgn_right,

        set_var,

        is_rational,

        refine_once,
        refine,

        other,             # PTOPO: characterize the param as double, cusp, pole etc
        get_sgn_deriv,

        display;


    local
        init,
        sgnL; 			   # Sign on left endpoint.

    option
        load = init;

        # Assumptions:
        #    1. We make the strong assumption that the
        #       defining polynomial is irreducible over the
        #       rationals. This is ok for maple since the function
        #       factors, provides such a functionallity.

        #
        # Initialization.
        #
        init := proc()
        local ta;
            poly := polynomial; # assure correct description.
            XX   := indets(poly)[1];
            sJ   := interval;
            J    := interval;
            mult := multiplicity;
            other := "";
            # print("here", poly, degree(poly));
            # if it is a rational number
            if (degree(poly, XX) = 1) then
                ta := solve(poly);
                sJ := [ta, ta];
                J := sJ;
            fi;
            if ( J[1] = J[2] ) then
                sgnL := 0;
            else
                # the isolating exI should not contain 0 inside
                sgnL := signum(eval(poly, XX = J[1]));
            fi;
            return;
        end proc; # init


        # #
        # # Methods
        # #
        get_degree := () -> degree(poly);
        get_interval := () -> J;
        get_left :=  ()-> J[1];
        get_right := ()-> J[2];
        get_approximation := ()-> evalf((J[1]+J[2])/2);
        get_midpoint := () -> (J[1]+J[2])/2;


        get_sgn_left  := () -> sgnL;
        get_sgn_right := () -> -sgnL;

        # # the sign of the derivative of the defining polynomial
        # # Notice that sign( deriv) = sign( right) - sign( left)
        # #                          = - sign( left)
        # getSignDeriv := () -> -sgnL;  # the sign of the derivative;

        set_var := proc(iXX)
            poly := subs(XX=iXX, poly);
            XX := iXX;
        end proc;

        is_rational := () -> `if`( J[1] = J[2], 'true', 'false');


        refine_once := proc()
        local m, fm;
            if is_rational() then return; fi;

            m := get_midpoint();
            fm := subs(XX = m, poly);
            if signum(fm) * get_sgn_left() < 0 then J := [J[1], m];
            else J := [m, J[2]]; fi;
            return;
        end proc;

        refine := proc(nb::posint := 1)
        local ii;
            for ii from 1 to nb do refine_once(); od;
        end proc;

        display := ()->
             printf( "< %a, [%a, %a], %5.10f >  mult=%a\n",
                     sort(poly, XX), J[1], J[2], evalf((J[1]+J[2])/2), mult);



            # the sign of the derivative of the defining polynomial
            # Notice that sign( deriv) = sign( right) - sign( left)
            #                          = - sign( left)
            get_sgn_deriv := () -> -sgnL;  # the sign of the derivative;


## Initialization
        init();
    end module;  # end rootof
end proc;        # end rootof
