function [x,v] = ploteigenvector (L,ev,ne,nsub,scale,fig)
%% declare local variables here if required by language
nv=ne*nsub + 1; 
Le=L/ne; 
dx=Le/nsub;
k=1;
x= zeros(nv);
v=zeros(nv); % declare and set to zero plot arrays
for e = 1:ne  % loop over elements
    xi= Le*(e-1); 
    vi= ev(2*e-1); 
    qi= ev(2*e); 
    vj= ev(2*e + 1); 
    qj= ev(2*e+2);
    for n= 1:nsub % loop over subdivisions
        xk=xi+dx*n; 
        x_scalar=(2*n-nsub)/nsub; % isoP coordinate
        vk=scale*(0.125*(4*(vi+vj)+2*(vi-vj)*(x_scalar^2-3)*x_scalar + Le*(x_scalar^2-1)*(qj-qi+(qi+qj)*x_scalar))); %Hermitian interpolant
        k = k+1;
        x(k) = xk;
        v(k)=vk; % build plot functions
    end  % end n loop
end % end e loop

end