
`evalr/pow` := proc(rg, fract, exact)
local ispos1, ispos2, newrg, newerrg, res, i, j, t1, t2,
    x,y,a,b,m;

    #print( "AAA"); print  ( type(rg,range) );
    if type(rg,range) then
        #print( "here", rg);

                
        if type(fract,integer) and 0 <= fract and irem(fract,2) <> 0 then
            return op(1,rg)^fract .. op(2,rg)^fract
        end if;
        ispos2 := `evalr/ispos`(op(2,rg),false,false,exact);
        if ispos2 = FAIL then
            return -infinity .. infinity
        elif ispos2 = false or ispos2 = 0 then
            if not type(fract,integer) then
                error "not a real number"
            else
                if exact then
                    return op(2,rg)^fract .. op(1,rg)^fract
                else
                    return evalf(op(2,rg)^fract) .. evalf(op(1,rg)^fract)
                end if
            end if
        elif ispos2 <> true then
            if op(1,rg) <> op(2,rg) then
                newrg := `evalr/chvar`(rg,op(1,ispos2));
                if not type(fract,integer) then
                    res := []
                else
                    t1 := `evalr/chvar`(rg,op(2,ispos2));
                    res := map(proc (x, fract) op(2,x)^fract .. op(1,x)^fract end proc,t1,fract)
                end if
            else
                if not type(fract,integer) then
                    t1 := `evalr/chvar`(op(1,rg),op(1,ispos2));
                    if exact then
                        t1^fract
                    else
                        evalf(t1^fract)
                    end if;
                    return t1 .. t1
                else
                    t1 := `evalr/chvar`(op(1,rg),ispos2);
                    return op(map(proc (x, fract, y) local t1; if y then t1 := x^fract else t1 := evalf(x^fract) end if; t1 .. t1 end proc,[op(op(1,t1)), op(op(2,t1))],fract,exact))
                end if
            end if
        else
            if op(1,rg) <> op(2,rg) then
                newrg := [rg];
                res := []
            else
                if exact then
                    return op(1,rg)^fract .. op(2,rg)^fract
                else
                    return evalf(op(1,rg)^fract) .. evalf(op(2,rg)^fract)
                end if
            end if
        end if;
        # print( "newrg", newrg);

        ## Power Hack 
        # print( "I am here");
        x := op(newrg): y := fract;
        if is( type( fract, posint)) then
            a := op(x)[1]; b := op(x)[2];
            # print( a, b);
            m := max( abs(a), b);
            if ( y mod 2 = 0 ) then
                res := [ op(res), a * m^(y-1) .. m^y];
            else                        
                res := [ op(res), a * m^(y-1) .. b * m^(y-1)];
            fi;
            if exact then
                return op(res)
            else
                return op(evalf(res))
            end if
        fi;
        
        for i to nops(newrg) do
            ispos1 := `evalr/ispos`(op([i, 1],newrg),false,false,exact);
            #print ("ispos1", ispos1);
            if ispos1 = FAIL then
                return -infinity .. infinity
            elif ispos1 = true then
                if 0 <= fract then
                    res := [op(res), op([i, 1],newrg)^fract .. op([i, 2],newrg)^fract]
                else
                    res := [op(res), op([i, 2],newrg)^fract .. op([i, 1],newrg)^fract]
                end if;
                next
            elif ispos1 = 0 then
                if 0 <= fract then
                    res := [op(res), 0 .. op([i, 2],newrg)^fract]
                else
                    res := [op(res), op([i, 2],newrg)^fract .. infinity]
                end if;
                next
            elif ispos1 <> false then
                if 0 <= fract then
                    res := [op(res), op(map(proc (x, fract) op(1,x)^fract .. op(2,x)^fract end proc,`evalr/chvar`(op(i,newrg),op(1,ispos1)),fract))]
                else
                    res := [op(res), op(map(proc (x, fract) op(2,x)^fract .. op(1,x)^fract end proc,`evalr/chvar`(op(i,newrg),op(1,ispos1)),fract))]
                end if;
                if type(fract,fraction) then
                    next
                end if;
                newerrg := `evalr/chvar`(op(i,newrg),op(2,ispos1))
            else
                if type(fract,fraction) then
                    if 0 <= fract then
                        res := [op(res), 0 .. op([i, 2],newrg)^fract]
                    else
                        res := [op(res), infinity .. op([i, 2],newrg)^fract]
                    end if;
                    next
                else
                    newerrg := [op(i,newrg)]
                end if
            end if;

            for j to nops(newerrg) do
                if fract < 0 and irem(fract,2) <> 0 then
                    res := [op(res), -infinity .. op([j, 1],newerrg)^fract, op([j, 2],newerrg)^fract .. infinity];
                    next
                end if;
                t1 := op(`evalr/evalr`(op([j, 1],newerrg)+op([j, 2],newerrg),exact));
                t2 := `evalr/ispos`(op(1,t1),false,false,exact);
                if t2 = FAIL then
                    return -infinity .. infinity
                elif t2 = true or t2 = 0 then
                    if 0 <= fract then
                        res := [op(res), 0 .. op([j, 2],newerrg)^fract]
                    else
                        res := [op(res), op([j, 2],newerrg)^fract .. infinity]
                    end if;
                    next
                end if;
                if op(1,t1) <> op(2,t1) then
                    t2 := `evalr/ispos`(op(2,t1),false,false,exact)
                end if;
                if t2 = false or t2 = 0 then
                    if 0 <= fract then
                        res := [op(res), 0 .. op([j, 1],newerrg)^fract]
                    else
                        res := [op(res), op([j, 1],newerrg)^fract .. infinity]
                    end if
                else
                    return -infinity .. infinity
                end if
            end do
        end do;
        if exact then
            op(res)
        else
            op(evalf(res))
        end if
    else
        t1 := map(`evalr/pow`,rg,fract,exact);
        if hastype(lt1,'infinity') then
            `simplify/infinity`(t1)
        else
            t1
        end if;
        if 1 < nops(t1) then
            `evalr/union`(t1)
        else
            t1
        end if
    end if
end proc;
