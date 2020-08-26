s = rng(100);
K = 536;             % Message length in bits, including CRC, K > 30
E = 1024;            % Rate matched output length, E <= 8192
EbNoArray = [0:0.2:3];         % EbNo in dB
Fer=zeros(1,length(EbNoArray)); %array of FER(snr)
L = 32;              % List length, a power of two, [1 2 4 8]
fe = 20;     % Number of frames to simulate
linkDir = 'UL';
if strcmpi(linkDir,'DL')
    % Downlink scenario (K >= 36, including CRC bits)
    crcLen = 24;      % Number of CRC bits for DL, Section 5.1, [6]
    poly = '24C';     % CRC polynomial
    nPC = 0;          % Number of parity check bits, Section 5.3.1.2, [6]
    nMax = 9;         % Maximum value of n, for 2^n, Section 7.3.3, [6]
    iIL = true;       % Interleave input, Section 5.3.1.1, [6]
    iBIL = false;     % Interleave coded bits, Section 5.4.1.3, [6]
else
    % Uplink scenario (K > 30, including CRC bits)
    crcLen = 11;
    poly = '11';
    nPC = 0;
    nMax = 10;
    iIL = false;
    iBIL = true;
end
R = K/E;                          % Effective code rate
bps = 2;                          % bits per symbol, 1 for BPSK, 2 for QPSK
for i=1:length(EbNoArray)
    EbNo=EbNoArray(i);
    EsNo = EbNo + 10*log10(bps);
    snrdB = EsNo + 10*log10(R);       % in dB
    noiseVar = 1./(10.^(snrdB/10));
    % Channel
    chan = comm.AWGNChannel('NoiseMethod','Variance','Variance',noiseVar);
    numferr = 0;
    numFrames=0;
    while numferr<fe
        numFrames=numFrames+1;
        % Generate a random message
        msg = randi([0 1],K-crcLen,1);
        % Attach CRC
        msgcrc = nrCRCEncode(msg,poly);
        % Polar encode
        encOut = nrPolarEncode(msgcrc,E,nMax,iIL);
        N = length(encOut);
        % Rate match
        modIn = nrRateMatchPolar(encOut,K,E,iBIL);
        % Modulate
        modOut = nrSymbolModulate(modIn,'QPSK');
        % Add White Gaussian noise
        rSig = chan(modOut);
        % Soft demodulate
        rxLLR = nrSymbolDemodulate(rSig,'QPSK',noiseVar);
        % Rate recover
        decIn = nrRateRecoverPolar(rxLLR,K,N,iBIL);
        % Polar decode
        decBits = nrPolarDecode(decIn,K,E,L,nMax,iIL,crcLen);
        % Compare msg and decoded bits
        numferr = numferr + any(decBits(1:K-crcLen)~=msg);

    end
    Fer(i)=numferr/numFrames;
    disp(['Block Error Rate: ' num2str(numferr/numFrames) ...
          ', at SNR = ' num2str(snrdB) ' dB'])
end
rng(s); 