clc;clear;

X= [1 3 5 ;2 4 6];
[M,N] = size(X);
% 法1：BasicDefinitionCumulant4的E1
result1 = zeros(M^2,M^2);
mu4 = zeros(M, M, M, M); % 用于存储所有(k1,k2,l1,l2)的结果
for k1 = 1:M
    for k2 = 1:M
        for l1 = 1:M
            for l2 = 1:M    
                
                temp = 0;
                for t = 1:N
                    temp = temp + X(k1, t) * conj(X(k2, t)) * X(l1, t) * conj(X(l2, t));
                end
                mu4(k1, k2, l1, l2) = temp / N;
                
                idx1 = (k1-1)*M + l1;
                idx2 = (k2-1)*M + l2;
                result1(idx1,idx2) = temp / N ;
                
            end
        end
    end
end
                              
% 法2：SingleSampleCumulant4的E1
result2 = zeros(M^2,M^2);
for i = 1: N
    xi = X(:,i);
    result2 = kron(xi,conj(xi)) * kron(xi,conj(xi))' / N + result2;
end
 
% 法3：holeSampleCumulant4的E1
XX = kron(X, conj(X));    
result3 = XX * XX';             
result3 = result3 / N;              
