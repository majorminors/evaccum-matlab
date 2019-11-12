function resp = idealObserver(intensity,curve)

         %get the estimated pcorrect
         id_p  =dsearchn(curve(:,1),intensity);
         
         p_int = curve(id_p,2);
         
         %random sample a response
         resp = randsample([0 1],1,true,[1-p_int p_int]);