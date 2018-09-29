
warning off
% load ('jaffe_213n_676d_10c_uni.mat');
% load ('jaffe_213n_676d_10c_uni_kernel_gaussian_0.01_post_Sample-Scale.mat');
% you should tune these parameters to obtain the best performance.
para1=[ 10 ];
para2=[ .0001   ];
para3=[  .001];


for i=1:length(para1)
    for j=1:length(para2)
        for k=1:length(para3)

            fprintf('params%12.6f%12.6f %12.6f%12.6f\n',para1(i),para2(j),para3(k))
            
        [result]=multiple3(K,y,para1(i),para2(j),para3(k))
        dlmwrite('yale.txt',[para1(i),para2(j),para3(k),result],'-append','delimiter','\t','newline','pc');

        end
    end
end
