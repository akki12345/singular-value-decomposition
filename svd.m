 m = input('Enter number of rows: '); 
 n = input('Enter number of columns: '); 
 
 for i = 1:m 
 for j = 1:n 
  str = ['Enter element in row ' num2str(i) ', col ' num2str(j) ': ']; 
  a(i,j) = input(str); 
 end 
 end 
 a
 B=size(a);   %In this code, we assume that a previous matrix "A" has already    been inputted.
 for j=1:B(1)
    for i=1:B(2)
        C(i,j)=a(j,i);
    end      %The transposed A-matrix should be C 
 end
 C
 
 [nRow1,nCol1]=size(a);
 [nRow2,nCol2]=size(C);
 
 if nCol1 ~= nRow2
 error('inner dimensions must match');
 end
 
 D = zeros(nRow1,nCol2);
 
 for i = 1:nRow1
    for j = 1:nCol2
        for k = 1:nCol1
 D(i,j) = D(i,j) + a(i,k)*C(k,j);
        end
    end
 end
D
d = D(1,1)*D(2,2) - D(1,2)*D(2,1);
d
t = D(1,1) + D(2,2);
e1 = (t + sqrt(t^2 - 4*d))/2;
e2 = (t - sqrt(t^2 - 4*d))/2;
if a(1,2) ~= 0
   x1 = [a(1,2); e1-a(1,1)];
   x2 = [a(1,2); e2-a(1,1)];
elseif a(2,1) ~= 0
   x1 = [e1-a(2,2); a(2,1)];
   x2 = [e2-a(2,2); a(2,1)];


else
   x1 = [1; 0];
   x2 = [0; 1];

end

disp(' ')
disp('For this matrix, the polynomial whose roots are the eigenvalues is:')
disp(['   e^2 - ' num2str(t) '*e + ' num2str(d) ' = 0'])
 
disp(' ')
disp('The first eigenvalue and eigenvector are:')
e1
x1
 
disp(' ')
disp('The second eigenvalue and eigenvector are:')
e2
x2
disp('')
disp('V of the matrix')
V=[x1,x2];
V
Z=size(V);   %In this code, we assume that a previous matrix "A" has already been inputted.
for j=1:Z(1)
    for i=1:Z(2)
        M(i,j)=V(j,i);
    end      %The transposed A-matrix should be C 
end
disp('')
disp('transpose of V')
M
disp('')
disp('singular values is:')
                                     
a1=sqrt(e1);
a2=sqrt(e2);
a1
a2
disp('')
disp('W of the matrix')
W=[a1,0;0,a2];
W
disp('')
disp('construct U of the matrix')
u1 = a*x1/a1;
u2 = a*x2/a2;
u1
u2;
U=[u1,u2];
U
disp('')
disp('finally SVD of the matrix a')
K = W*M;
  [nRow1,nCol1]=size(U);
  [nRow2,nCol2]=size(K);
 
 if nCol1 ~= nRow2
 error('inner dimensions must match');
 end
Q = zeros(nRow1,nCol2);

for i = 1:nRow1
    for j = 1:nCol2
        for k = 1:nCol1
 Q(i,j) = Q(i,j)+ U(i,k)*K(k,j);
        end
    end
 end
 Q   

