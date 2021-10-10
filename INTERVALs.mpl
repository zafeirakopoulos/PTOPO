INTERVALs := module()
export `+`, `-`, `*`, `/`, `^`, eval,
    ln, exp,
    sin, cos, tan, sec, csc, cot,
    arcsin, arccos, arctan, arcsec, arccsc, arccot,
    sinh, cosh, tanh, sech, csch, coth,
    arcsinh, arccosh, arctanh, arcsech, arccsch, arccoth,
    signum, NewProd;
option package;
local combine, i;

    `+` := proc(x, y)
        if nargs = 0 then
            0;
        elif nargs = 1 then
            x;
        elif nargs = 2 then
            if type( x, 'specfunc'('`..`', 'INTERVAL') ) or
            type( y, 'specfunc'('`..`', 'INTERVAL') ) then
                evalr( x + y );
            else
                x + y;
            end if;
        else
            procname( procname( x, y ), args[3..-1] );
        end if;
    end proc;

    `-` := proc(x, y)
        if nargs = 1 then
            if type( x, 'specfunc'('`..`', 'INTERVAL') ) then
                INTERVAL( seq( -op( 2, i )..-op( 1, i ), i = x ) );
            else
                -x;
            end if;
        elif nargs = 2 then
            `+`( x, `-`(y) );
        else
            error "expecting 1 or 2 arguments, but got %1", nargs;
        end if;
    end proc;


    NewProd := proc( x, y)
    local a, b, c, d, t;
        a := op( op( x))[1]; b := op( op( x))[2];
        c := op( op( x))[1]; d := op( op( x))[2];
        
        t := a*c, a*d, b*c, b*d;
        return INTERVAL( min( t).. max( t));
    end proc;
    
   `*` := proc(x, y)
        if nargs = 0 then
            1;
        elif nargs = 1 then
            x;
        elif nargs = 2 then
            if type( x, 'specfunc'('`..`', 'INTERVAL') ) or
            type( y, 'specfunc'('`..`', 'INTERVAL') ) then
               NewProd( x, y); 
               #evalr( x * y );
            else
                x * y;
            end if;
        else
            procname( procname( x, y ), args[3..-1] );
        end if;
    end proc;

    `/` := proc(x, y)
        if nargs = 1 then
            if type( x, 'specfunc'('`..`', 'INTERVAL') ) then
                evalr( 1/x );
            else
                1/x;
            end if;
        elif nargs = 2 then
            `*`( x, `/`(y) );
        else
            error "expecting 1 or 2 arguments, but got %1", nargs;
        end if;
    end proc;

    `^` := proc(x, y)
        local a, b, m;
        if nargs = 2 then
            if type( x, 'specfunc'('`..`', 'INTERVAL') ) or
            type( y, 'specfunc'('`..`', 'INTERVAL') ) then
               if is( type( y, posint)) then
                    a := op( op( x))[1]; b := op( op( x))[2]; 
                    m := max( abs(a), b);
                    if ( y mod 2 = 0 ) then
                        INTERVAL( a * m^(y-1) .. m^y);
                    else                        
                        INTERVAL( a * m^(y-1) .. b * m^(y-1));
                    fi;
                else
                    evalr( x ^ y );
                fi;
                #evalr( x ^ y );
            else
                x ^ y;
            end if;
        else
            error "expecting 2 arguments, but got %1", nargs;
        end if;
    end proc;

    eval := proc( x::uneval, y )
        if nargs = 1 then
            :-eval( ':-eval'( x ) );
        elif type( y, 'posint' ) then
            :-eval( ':-eval'( x, y ) );
        elif type( y, {'`=`', 'list'('`=`')} ) then
            try
                combine( :-eval(x), y );
            catch:
                error;
            end try;
        else
            error "wrong number (or type) of parameters in function eval";
        end if;
    end proc;

    combine := proc(x, y)
    local M, i, result;

        if type( x, atomic ) then
            :-eval( x, y );
        elif type( x, 'specfunc'('`..`', INTERVAL) ) then
            :-eval( x, y );
        elif type( x, ':-`*`' ) then
            result := seq( procname( i, args[2..-1] ), i = x );

            `*`( result );
        elif type( x, ':-`+`' ) then
            result := seq( procname( i, args[2..-1] ), i = x );

            `+`( result );
        elif type( x, ':-`^`' ) then
            result := seq( procname( i, args[2..-1] ), i = x );

            `^`( result );
        elif type( x, 'function' ) then
            evalr( :-eval( op( 0, x ), y )( seq( procname( i, args[2..-1] ), i = x ) ) );
        else
            map( procname, x, args[2..-1] );
        end if;
    end proc;

    for i in [
        ln, exp,
        sin, cos, tan, sec, csc, cot,
        arcsin, arccos, arctan, arcsec, arccsc, arccot,
        sinh, cosh, tanh, sech, csch, coth,
        arcsinh, arccosh, arctanh, arcsech, arccsch, arccoth
             ] do
        assign( i = FromInert(
            _Inert_PROC(
                _Inert_PARAMSEQ( _Inert_NAME("x") ),
                _Inert_LOCALSEQ(),
                _Inert_OPTIONSEQ(),
                _Inert_EXPSEQ(),
                _Inert_STATSEQ(
                    _Inert_FUNCTION(
                        _Inert_NAME("evalr"),
                        _Inert_EXPSEQ(
                            _Inert_FUNCTION(
                                _Inert_MEMBER(
                                    _Inert_EXPSEQ(),
                                    _Inert_NAME(convert(i, 'string'))
                                             ),
                                _Inert_EXPSEQ(
                                    _Inert_PARAM(1)
                                             )
                                           )
                                     )
                                   )
                              ),
                _Inert_DESCRIPTIONSEQ(),
                _Inert_GLOBALSEQ(),
                _Inert_LEXICALSEQ()
                       )
                             ) );
    end do;

    signum := proc(x)
    local result;

        if nargs = 1 then
            result := evalr( Signum( x ) );

            if result = FAIL then
                undefined;
            else
                result;
            end if;
        else
            error "not yet implemented";
        end if;
    end proc;
end module:

with(INTERVALs):
