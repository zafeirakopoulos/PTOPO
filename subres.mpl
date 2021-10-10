
with(LinearAlgebra):



PRS := module()

export
    Var,
    StHa,
    display;



local
    poly_lc_sgn;


    poly_lc_sgn := proc ( f )
        signum( lcoeff( f ) );
    end proc:



#############################################
# Find Variance on an evaluation (at a point) #
# of a Polynomial Remainder Sequence.         #
#############################################


    Var := proc (L)
        # Local variables.
    local i, temp, SIGN, variations;

        # Check for empty list.
        if ( nops(L) = 0 ) then
            return 0
        end if;

        #
        # Procceed for computation.
        #
        temp := op(1, L);
        if ( poly_lc_sgn ( temp ) > 0 ) then
            SIGN := 1
        elif ( poly_lc_sgn ( temp ) = 0 ) then
            SIGN := 0
        else
            SIGN := -1
        end if;
        variations := 0;
        for i from 2 to nops(L) do
            temp := op(i, L);
            if ( SIGN = 0 ) and ( poly_lc_sgn ( temp ) > 0 ) then
                SIGN := 1
            elif ( SIGN = 0 ) and ( poly_lc_sgn ( temp ) < 0 ) then
                SIGN := -1
            end if;
            if ( poly_lc_sgn ( temp ) > 0 ) and ( SIGN < 0 ) then
                variations := variations + 1;
                SIGN := 1
            elif ( poly_lc_sgn ( temp ) < 0 ) and ( SIGN > 0 ) then
                variations := variations + 1;
                SIGN := -1
            end if;
        end do;
        variations;
    end proc:




    # compute the Sylvester-Habicht sequence
    StHa:= proc(a, b, x)
    local suc,p,q,n1,n2,n,m,lp,lq,r,R,S0,S1,d,NS,Sr,j;
        #p:=expand(a); q:=expand(diff(a,x)*b);
        p := a;
        q := b;
        if q=0 then return x,p fi;
        n1:=degree(p,x); n2:=degree(q,x); lp:=lcoeff(p,x); lq:=lcoeff(q,x);
        suc:=p,q;

        if n1>n2 then
            n:=n1-1; d:=n-n2;
            if d<>0 then
                S1:=((-1)**(d*(d+1)/2))*(lp*lq)**d*q; S0:=((-1)**((d+2)*(d+3)/2))*lp**d*prem(p,q,x);
                if S0<>0 then suc:=suc,sort( S1, x), sort( S0, x) else suc:=suc,sort( S1, x) fi;
            else S1:=q; S0:= -prem(p,q,x); if S0<>0 then suc:= suc, sort( S0, x) fi;
            fi;
            j:=n2-1;
        else
            n:=n2; S1:=q; S0:= -rem(lq**2*p,q,x); j:=n-1; if S0<>0 then suc:=suc,sort( S0,x) fi;
        fi;
        R:=lcoeff(S1,x); if S0=0 then r:=-1 else r:=degree(S0,x) fi;
        while r>=0 do
            if r>0 then
                if r=j then
                    NS:='NS'; divide(prem(S1,S0,x),((-R)**(j-r+2)),NS); NS:=expand(NS)*((-1)**((j-r+2)*(j-r+3)/2));
                    if NS<>0 then suc:=suc,sort( NS, x) fi;
                    S1:=S0; S0:=NS; R:=lcoeff(S1,x); j:=r-1;
                else
                    for m from r+1 to j-1 do suc:=suc,0 od: Sr:='Sr';
                    divide(prem(S1,S0,x),(R**(j-r+2)),Sr); S1:='S1';
                    divide(((lcoeff(S0,x)**(j-r))*S0),(R**(j-r)),S1);
                    S1:=S1*((-1)**((j-r)*(j-r+1)/2)); S0:=Sr*((-1)**((j-r+2)*(j-r+3)/2));
                    if S0<>0 then suc:=suc,sort(S1, x), sort( S0, x) else suc:=suc, sort( S1, x) fi;
                    R:=lcoeff(S1,x); j:= r-1;
                fi
            else
                if j>0 then
                    S1:='S1'; divide(((lcoeff(S0,x)**j)*S0),(R**j),S1);
                    S1:=S1*((-1)**(j*(j+1)/2)); suc:=suc,sort( S1, x);
                fi:
                S0:=0;
            fi:
            if S0<>0 then r:=degree(S0,x) else r:=-1 fi:
        od:
        [ suc ] ;
    end proc;


    display := proc(L)
    local i, dd;
        dd := nops(L) - 1;
        for i from 1 to nops( L) do
            printf( "S[%2d] : %a \n", dd - i + 1, L[i]);
            # printf( "S[%2d] : %a \n", dd - i + 1, L[i]);
        od;
    end proc;



end module: # end PRS
