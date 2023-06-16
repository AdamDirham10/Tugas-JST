function output = LMFnlsq( fungsi, data, epsX, epsFX, jacobian, maksIterasi, minBobotSelisih, maksBobotSelisih)

jumlahData = length(data);
epsX       = epsX * ones(jumlahData,1);

lambda     = 0;
lambdaC    = 1;

outputData = feval(fungsi,data);

jumlahKuadratOutputData = outputData'*outputData;

[JTJ,JTOutput] = HitungJTJJTOutput(fungsi,jacobian,data,outputData,epsX);

skalaD = diag(JTJ);
skalaD(skalaD<=0)=1;
akarSkalaD = sqrt(skalaD);

iterasi = 0;
iterasiDeltaJumlahKuadrat = 1;

while 1
    iterasi = iterasi+1;
    deltaS = zeros(jumlahData,1);
    
    diagJTJ = diag(JTJ);
    
    while 1

        while 1
            UA = triu(JTJ,1);
            JTJ = UA' + UA + diag(diagJTJ+lambda*skalaD);
            [U,p] = chol(JTJ);            
            if p==0, break, end
            
            lambda = 2*lambda;
            if lambda==0, lambda=1; end
        end

        deltaX = U\(U'\JTOutput);
        bobotJTOutput = deltaX'*JTOutput;
        if bobotJTOutput<=0, break, end

        for i=1:jumlahData
            z = diagJTJ(i)*deltaX(i);
            if i>1, z=JTJ(i,1:i-1)*deltaX(1:i-1)+z; end
            if i<jumlahData, z=JTJ(i+1:jumlahData,i)'*deltaX(i+1:jumlahData)+z; end            
            deltaS(i) = 2*JTOutput(i)-z;
        end
        
        deltaQ = deltaS'*deltaX;
        
        deltaS = data-deltaX;
        outputTerhitung = feval(fungsi,deltaS);

        iterasiDeltaJumlahKuadrat = iterasiDeltaJumlahKuadrat+1;
        
        jumlahKuadratOutputTerhitung = outputTerhitung'*outputTerhitung;
        deltaOutput  = jumlahKuadratOutputData-jumlahKuadratOutputTerhitung;
        
        isSelesai = 1;
        if all((abs(deltaX)-epsX)<=0) || iterasiDeltaJumlahKuadrat>=maksIterasi || abs(deltaOutput)<=epsFX
            break
        end        
        isSelesai = 0;
        
        if deltaOutput>=minBobotSelisih*deltaQ, break, end
        
        JTJ = U;
        y = .5;
        z = 2*bobotJTOutput-deltaOutput;
        if z>0, y=bobotJTOutput/z; end
        if y>.5, y=.5; end
        if y<.1, y=.1; end
        
        if lambda==0

            y = 2*y;
            
            for i = 1:jumlahData
                JTJ(i,i) = 1/JTJ(i,i);
            end
            
            for i = 2:jumlahData
                ii = i-1;
                for j= 1:ii
                    JTJ(j,i) = -JTJ(j,j:ii)*JTJ(j:ii,i).*JTJ(i,i);
                end
            end
            
            for i = 1:jumlahData
                for j= i:jumlahData
                    JTJ(i,j) = abs(JTJ(i,j:jumlahData)*JTJ(j,j:jumlahData)');
                end
            end
            
            lambda = 0;
            for i = 1:jumlahData
                z = JTJ(1:i,i)'*akarSkalaD(1:i)+z;
                if i<jumlahData
                    ii = i+1;
                    z  = JTJ(i,ii:jumlahData)*akarSkalaD(ii:jumlahData)+z;
                end
                z = z*akarSkalaD(i);
                if z>lambda, lambda=z; end
            end
            
            tr = diag(JTJ)'*skalaD;
            if tr<lambda, lambda=tr; end
            lambda  = 1/lambda;
            lambdaC = lambda;
        end
        
        lambda = lambda/y;
        
        if deltaOutput>0, deltaOutput=-1e300; break, end
    end
    
    if isSelesai, break, end
    
    if deltaOutput > maksBobotSelisih*deltaQ
        lambda = lambda/2;
        if lambda < lambdaC, lambda=0; end
    end
    
    data=deltaS;
    outputData=outputTerhitung;
    jumlahKuadratOutputData=jumlahKuadratOutputTerhitung;
    [JTJ,JTOutput] = HitungJTJJTOutput(fungsi,jacobian,data,outputData,epsX);    
end

output = data;

function [JTJ,JTOutput] = HitungJTJJTOutput(fungsi,jacobian,data,outputData,epsX)
if isa(jacobian,'function_handle')
    J = jacobian(data);
else
    J = feval(jacobian,fungsi,outputData,data,epsX);
end
JTJ = J'*J;
JTOutput = J'*outputData;
 
function J = finiteJacobian(fungsi,outputData,data,epsX)
jumlahData = length(data);
J  = zeros(length(outputData),jumlahData);

for k = 1:jumlahData
    deltaX = .25 * epsX(k);
    if deltaX==0, deltaX=eps; end
    xData = data;
    xData(k) = xData(k) + deltaX;
    outputTerhitung = feval(fungsi,xData);
    
    J(:,k)=((outputTerhitung(:)-outputData(:))/deltaX);
end