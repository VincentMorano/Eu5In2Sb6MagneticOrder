function [f2,f2Err] = structFact(R0,M,int,intErr,q)
%structFact Calculate structure factors.
%   Given R0 and M from ResLib, the integrated intensity, its error, and q,
%   return the structure factor for an a3 scan.

M22=M(2,2,:);
M22=reshape(M22,[length(M22),1,1]);
f2=q.*int.*sqrt(M22./(2*pi))./R0';
f2Err=q.*intErr.*sqrt(M22./(2*pi))./R0'; % See Jonathan's Mn3Sn supplement
end