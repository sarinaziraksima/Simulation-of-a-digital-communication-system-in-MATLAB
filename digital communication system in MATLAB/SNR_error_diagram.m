


clear;
%close all;

%%%%%%%%%%%%%%% note:
% r is sent data bits
% r1 is domain of sent data
% a is recived data bits
% a1 is domain of recived data

%%%%%%%%%%%%%%% initial values
M=2;  % M = Number of layers example: m=2 -> [-1 1] m=3 -> [-2 0 2] m=4 -> [-3 -1 1 3] 
n=1000; % n = Number of bits generated
Tb=1; % Tb = Duration of each bit sent by a pulse
delta = 0.01;
snr= -10:2:10;
error = zeros(1,length(snr));
%%%%%%%%%%%%%%% Channal initial values

% if channel is :   alfa.exp(-(1/beta)tx);
alfa = 0.1; % domain of channel
beta = 0.03; % frequency of channel is 1/beta = 1k
gain = 1;



for l=1:length(snr)
    

    %%%%%%%%%%%%%%% sending data generation


    
rr = randi([0 1],1,n);


%%%%%%%%%%%%%%% channel coding

k = 4;
nn=n+mod(n,k);
rr(length(rr)+1:nn)=0;

P = [[1 0 1];[1 1 1];[0 1 1];[1 1 0]];
G = [eye(k) P];

HT = [P;eye(3)];


for i=1:length(rr)/k

    d(i,:)=rr((i-1)*k+1:i*k);
    
end
c=d*G;

c=mod(c,2);

%%%%%%%%%%%%%%% amplitute of sending data generation



for i=1:length(c(:,1))

    r((i-1)*7+1:i*7)=c(i,1:7);
    
end
r1=r;
r1(r1==0)=-1;
%%%%%%%%%%%%%%% analog pulse generation

tx = 0:delta:Tb*length(r1); 
y = ones(1,(1/delta)*Tb*length(r1));

t_ref = 0:delta:Tb;
 
Y_ref = [zeros(1,Tb*0.4/delta) ones(1,length(t_ref)-Tb*0.8/delta) zeros(1,Tb*0.4/delta)]; % Regtangular with duty cycle 20%
%Y_ref = ones(1,length(t_ref)); % Square pulse
%Y_ref = [zeros(1,40) ones(1,length(t_ref)-80) zeros(1,40)]; % Regtangular with duty cycle (Tb*100)-80/(Tb*100)-> if Tb=1s thus: duty cycle=0.2 or 20%
%Y_ref = [zeros(1,20) ones(1,length(t_ref)-40) zeros(1,20)]; % Regtangular with duty cycle (Tb*100)-40/(Tb*100)-> if Tb=1s thus: duty cycle=0.2 or 60%
%Y_ref = 2*[0:delta:Tb/2 Tb/2-delta:-delta:0]/Tb; % triangular pulse
%Y_ref = cos(t_ref*2*pi/Tb); % Sinusoidal pulse
%Y_ref = sawtooth(2*pi*t_ref*Tb); % Sawtooth pulse

for i = 1:1:length(r1)
    
    y(1,(Tb*(i-1)/delta)+1:Tb*i/delta+1) = Y_ref*r1(1,i);
    
end






%%%%%%%%%%%%%%% HT filter

    if 1==1
        
        h_c = [alfa*exp(-1*(1/beta)*tx) zeros(1,2*length(y))];
        y=[y zeros(1,2*length(y))];
        
        H_c = fft(h_c);
        
        HT=(sqrt(abs(H_c).^-1));
        clear g;
        g=gain*(ifft(HT.*fft(y)))*delta;
        
        g=g(1:length(tx));
       
        
        
        
        
       
       y=g;
        
        
    end
    
    



%%%%%%%%%%%%%%% Channel

    if 1==1
        Fs = 1/delta;
        
        h_c = [alfa*exp(-1*(1/beta)*tx) zeros(1,2*length(y))];
        y=[y zeros(1,2*length(y))];
        
        
        clear g;
        g=(ifft(fft(h_c).*fft(y)))*delta;
        
        g=g(1:length(tx));
       
        
        y=g;
    end
    
  
    


%%%%%%%%%%%%%%% adding noise
if 1==1
    y = awgn(y,snr(l),'measured');

    
end

%%%%%%%%%%%%%%% HR filter

    if 1==1
        
        
        final_amp = 5000000/4.954;
        Fs = 1/delta;
        
        h_c = [alfa*exp(-1*(1/beta)*tx) zeros(1,2*length(y))];
        y=[y zeros(1,2*length(y))];
        
        H_c = fft(h_c);
        
        HR=(sqrt(abs(H_c).^-1));
        clear g;
        g=final_amp*(ifft(HR.*fft(y)))*delta/gain;
        
        g=g(1:length(tx));
       
        
        y=g;
    end
    
    







%%%%%%%%%%%%%%% sampeling and digital resived data

for i = 1:1:length(r1)

    a1(1,i) = y(1,(i-0.5)*Tb/delta+1); 

end


if mod(M,2)==1
   a1 = round((a1)/2)*2; 

end
if mod(M,2)==0
    a1 = round((a1)/2+0.5)*2-1;
end

a1(a1>max(r1))=max(r1);
a1(a1<min(r1))=min(r1);

a=a1;
a(a==-1)=0;





for i=1:length(a)/7

    resived_c1(i,1:7)=a((i-1)*7+1:i*7);
    
end


%%%%%%%%%%%%%%% channel decoding
     
if 1==1
    k = 4;

    P = [[1 0 1];[1 1 1];[0 1 1];[1 1 0]];
    G = [eye(k) P];

    HT = [P;eye(3)];

    clear s; 
    s=resived_c1*HT;
    s=mod(s,2);

    for i=1:length(resived_c1(:,1)/length(rr))

       if(s(i)==[1 0 1])
           resived_c1(i,:)=resived_c1(i,:)+[1 0 0 0 0 0 0];
       end
       if(s(i)==[1 1 1])
           resived_c1(i,:)=resived_c1(i,:)+[0 1 0 0 0 0 0];
       end
       if(s(i)==[0 1 1])
           resived_c1(i,:)=resived_c1(i,:)+[0 0 1 0 0 0 0];
       end
       if(s(i)==[1 1 0])
           resived_c1(i,:)=resived_c1(i,:)+[0 0 0 1 0 0 0];
       end
       if(s(i)==[1 0 0])
           resived_c1(i,:)=resived_c1(i,:)+[0 0 0 0 1 0 0];
       end
       if(s(i)==[0 1 0])
           resived_c1(i,:)=resived_c1(i,:)+[0 0 0 0 0 1 0];
       end
       if(s(i)==[0 0 1])
           resived_c1(i,:)=resived_c1(i,:)+[0 0 0 0 0 0 1];
       end


    end

    resived_c1 = mod(resived_c1,2);
end
a=ones(1,length(rr));
for i=1:length(resived_c1(:,1))

    a((i-1)*4+1:i*4)=resived_c1(i,1:4);
    
end


rr=rr(1:n);
a=a(1:n);



    error(1,l) = sum(abs(rr-a)>0);



end

figure;
plot(snr,error*100/n)
title('Plot for Error')
ylabel('error%')
xlabel('snr')
