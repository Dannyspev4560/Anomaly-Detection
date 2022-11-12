function threshold = delta
    threshold=0;
    if delta< exp
        threshold=sqrt(2/pi)*(1/delta);
    else
        threshold=sqrt(2/(pi*exp(1)))*(1/log(delta)));
    end
endn=768;
delta = (287-1)/92;