function T=path_transitivity(W)
%TRANSITIVITY_WU    Transitivity
%
%   T = path_transitivity_wu(W);
%
%   Transitivity is the ratio of 'triangles to triplets' in the network.
%   (A classical version of the clustering coefficient).
%
%   Input:      W       weighted undirected connection matrix
%
%   Output:     T       transitivity scalar
%
%   Reference: Rubinov M, Sporns O (2010) NeuroImage 52:1059-69
%              based on Onnela et al. (2005) Phys Rev E 71:065103
%
%
%   Mika Rubinov, UNSW, 2010

%% shortest path
n=length(W);
%w=zeros(1,n);
m=zeros(n,n);
T=zeros(n,n);
dis=zeros(n,n);

%%
%w=strengths_und (W);
dis=1./W;

%% ---matching index---%% 

for i=1:n-1
    for j=i+1:n
        x=0;
        y=0;
        z=0;
       for k=1:n
           if W(i,k)~=0 && W(j,k)~=0 && k~=i && k~=j
               x=x+W(i,k)+W(j,k);
           end
           if k~=j
               y=y+W(i,k);
           end
           if k~=i
               z=z+W(j,k);
           end
       end
       m(i,j)=x/(y+z);
    end
end
   m=m+m';
 M=m;  
[D P B] = f_FastFloyd_und(dis);

%% --- path transitivity ---%% 
for i=1:n-1
    for j=i+1:n
        x=0;
        path = retrieve_shortest_path(i,j,P,B)
         K=length(path);
         
         for t=1:K-1
             for l=t+1:K
                 x=x+m(path(t),path(l));
             end
         end
        T(i,j)=2*x/(K*(K-1));      
    end
end
T=T+T';
