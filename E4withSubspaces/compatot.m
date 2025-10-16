load MA1PON

T1 = sqrt(mean((T11(7:9,1:3)').^2)');
T1 = [T1 sqrt(mean((T11(7:9,4:6)').^2)')];
T1 = [T1 sqrt(mean((T12(7:9,1:3)').^2)')];
T1 = [T1 sqrt(mean((T12(7:9,4:6)').^2)')];
T1 = [T1 sqrt(mean((T13(7:9,1:3)').^2)')];
T1 = [T1 sqrt(mean((T13(7:9,4:6)').^2)')];

load AR2PON

T1 = [T1 sqrt(mean((T11(7:9,1:5)').^2)')];
T1 = [T1 sqrt(mean((T11(7:9,6:10)').^2)')];
T1 = [T1 sqrt(mean((T12(7:9,1:5)').^2)')];
T1 = [T1 sqrt(mean((T12(7:9,6:10)').^2)')];

load MAX1PON

T1 = [T1 sqrt(mean((T11(7:9,1:5)').^2)')];
T1 = [T1 sqrt(mean((T11(7:9,6:10)').^2)')];
T1 = [T1 sqrt(mean((T12(7:9,1:5)').^2)')];
T1 = [T1 sqrt(mean((T12(7:9,6:10)').^2)')];

load ARX2PON

T1 = [T1 sqrt(mean((T11(7:9,1:8)').^2)')];
T1 = [T1 sqrt(mean((T11(7:9,9:16)').^2)')];
T1 = [T1 sqrt(mean((T12(7:9,1:8)').^2)')];
T1 = [T1 sqrt(mean((T12(7:9,9:16)').^2)')];

load MA1EXT

T2 = sqrt(mean((T11(7:9,1:3)').^2)');
T2 = [T2 sqrt(mean((T11(7:9,4:6)').^2)')];
T2 = [T2 sqrt(mean((T12(7:9,1:3)').^2)')];
T2 = [T2 sqrt(mean((T12(7:9,4:6)').^2)')];
T2 = [T2 sqrt(mean((T13(7:9,1:3)').^2)')];
T2 = [T2 sqrt(mean((T13(7:9,4:6)').^2)')];

load AR2EXT

T2 = [T2 sqrt(mean((T11(7:9,1:5)').^2)')];
T2 = [T2 sqrt(mean((T11(7:9,6:10)').^2)')];
T2 = [T2 sqrt(mean((T12(7:9,1:5)').^2)')];
T2 = [T2 sqrt(mean((T12(7:9,6:10)').^2)')];

load MAX1EXT

T2 = [T2 sqrt(mean((T11(7:9,1:5)').^2)')];
T2 = [T2 sqrt(mean((T11(7:9,6:10)').^2)')];
T2 = [T2 sqrt(mean((T12(7:9,1:5)').^2)')];
T2 = [T2 sqrt(mean((T12(7:9,6:10)').^2)')];

load ARX2EXT

T2 = [T2 sqrt(mean((T11(7:9,1:8)').^2)')];
T2 = [T2 sqrt(mean((T11(7:9,9:16)').^2)')];
T2 = [T2 sqrt(mean((T12(7:9,1:8)').^2)')];
T2 = [T2 sqrt(mean((T12(7:9,9:16)').^2)')];

load MA1OBS

T3 = sqrt(mean((T11(7:9,1:3)').^2)');
T3 = [T3 sqrt(mean((T11(7:9,4:6)').^2)')];
T3 = [T3 sqrt(mean((T12(7:9,1:3)').^2)')];
T3 = [T3 sqrt(mean((T12(7:9,4:6)').^2)')];
T3 = [T3 sqrt(mean((T13(7:9,1:3)').^2)')];
T3 = [T3 sqrt(mean((T13(7:9,4:6)').^2)')];

load AR2OBS

T3 = [T3 sqrt(mean((T11(7:9,1:5)').^2)')];
T3 = [T3 sqrt(mean((T11(7:9,6:10)').^2)')];
T3 = [T3 sqrt(mean((T12(7:9,1:5)').^2)')];
T3 = [T3 sqrt(mean((T12(7:9,6:10)').^2)')];

load MAX1OBS

T3 = [T3 sqrt(mean((T11(7:9,1:5)').^2)')];
T3 = [T3 sqrt(mean((T11(7:9,6:10)').^2)')];
T3 = [T3 sqrt(mean((T12(7:9,1:5)').^2)')];
T3 = [T3 sqrt(mean((T12(7:9,6:10)').^2)')];

load ARX2OBS

T3 = [T3 sqrt(mean((T11(7:9,1:8)').^2)')];
T3 = [T3 sqrt(mean((T11(7:9,9:16)').^2)')];
T3 = [T3 sqrt(mean((T12(7:9,1:8)').^2)')];
T3 = [T3 sqrt(mean((T12(7:9,9:16)').^2)')];

load MA1MET

T4 = sqrt(mean((T11(7:9,1:3)').^2)');
T4 = [T4 sqrt(mean((T11(7:9,4:6)').^2)')];
T4 = [T4 sqrt(mean((T12(7:9,1:3)').^2)')];
T4 = [T4 sqrt(mean((T12(7:9,4:6)').^2)')];
T4 = [T4 sqrt(mean((T13(7:9,1:3)').^2)')];
T4 = [T4 sqrt(mean((T13(7:9,4:6)').^2)')];

load AR2MET

T4 = [T4 sqrt(mean((T11(7:9,1:5)').^2)')];
T4 = [T4 sqrt(mean((T11(7:9,6:10)').^2)')];
T4 = [T4 sqrt(mean((T12(7:9,1:5)').^2)')];
T4 = [T4 sqrt(mean((T12(7:9,6:10)').^2)')];

load MAX1MET

T4 = [T4 sqrt(mean((T11(7:9,1:5)').^2)')];
T4 = [T4 sqrt(mean((T11(7:9,6:10)').^2)')];
T4 = [T4 sqrt(mean((T12(7:9,1:5)').^2)')];
T4 = [T4 sqrt(mean((T12(7:9,6:10)').^2)')];

load ARX2MET

T4 = [T4 sqrt(mean((T11(7:9,1:8)').^2)')];
T4 = [T4 sqrt(mean((T11(7:9,9:16)').^2)')];
T4 = [T4 sqrt(mean((T12(7:9,1:8)').^2)')];
T4 = [T4 sqrt(mean((T12(7:9,9:16)').^2)')];

load MA1MET

T5 = sqrt(mean((T11(7:9,2:3)').^2)');
T5 = [T5 sqrt(mean((T11(7:9,5:6)').^2)')];
T5 = [T5 sqrt(mean((T12(7:9,2:3)').^2)')];
T5 = [T5 sqrt(mean((T12(7:9,5:6)').^2)')];
T5 = [T5 sqrt(mean((T13(7:9,2:3)').^2)')];
T5 = [T5 sqrt(mean((T13(7:9,5:6)').^2)')];

load AR2MET

T5 = [T5 sqrt(mean((T11(7:9,3:5)').^2)')];
T5 = [T5 sqrt(mean((T11(7:9,8:10)').^2)')];
T5 = [T5 sqrt(mean((T12(7:9,3:5)').^2)')];
T5 = [T5 sqrt(mean((T12(7:9,8:10)').^2)')];

load MAX1MET

T5 = [T5 sqrt(mean((T11(7:9,4:5)').^2)')];
T5 = [T5 sqrt(mean((T11(7:9,9:10)').^2)')];
T5 = [T5 sqrt(mean((T12(7:9,4:5)').^2)')];
T5 = [T5 sqrt(mean((T12(7:9,9:10)').^2)')];

load ARX2MET

T5 = [T5 sqrt(mean((T11(7:9,6:8)').^2)')];
T5 = [T5 sqrt(mean((T11(7:9,14:16)').^2)')];
T5 = [T5 sqrt(mean((T12(7:9,6:8)').^2)')];
T5 = [T5 sqrt(mean((T12(7:9,14:16)').^2)')];
