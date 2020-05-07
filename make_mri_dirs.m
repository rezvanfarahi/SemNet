out_path = '/imaging/rf02/Semnet';

% sub-directories for individual subjects
subs_dir = {'/MRI_meg16_0030/mri/orig/', ... %0
            '/MRI_meg16_0032/mri/orig/', ... %1
            '/MRI_meg16_0033/mri/orig/', ... %2
            '/MRI_meg16_0034/mri/orig/', ... %3
            '/MRI_meg16_0035/mri/orig/', ... %4
            '/MRI_meg16_0039/mri/orig/', ... %5
            '/MRI_meg16_0041/mri/orig/', ... %6
            '/MRI_meg16_0042/mri/orig/', ... %7
            '/MRI_meg16_0045/mri/orig/', ... %8
            '/MRI_meg16_0047/mri/orig/', ... %9
            '/MRI_meg16_0052/mri/orig/', ... %10
            '/MRI_meg16_0056/mri/orig/',... %11
            '/MRI_meg16_0069/mri/orig/',... %12 
            '/MRI_meg16_0070/mri/orig/', ... %13
            '/MRI_meg16_0071/mri/orig/', ... %14
            '/MRI_meg16_0072/mri/orig/', ... %15
            '/MRI_meg16_0073/mri/orig/', ... %16
            '/MRI_meg16_0075/mri/orig/', ... %17
            '/MRI_meg16_0078/mri/orig/', ... %18
            '/MRI_meg16_0082/mri/orig/', ... %19
            '/MRI_meg16_0086/mri/orig/', ... %20
            '/MRI_meg16_0097/mri/orig/', ... %21 
            '/MRI_meg16_0122/mri/orig/', ... %22
            '/MRI_meg16_0123/mri/orig/', ... %23
            '/MRI_meg16_0125/mri/orig/', ... %24
                    
                   
    };


for cnt = 1:length(subs_dir)
    this_dir=[out_path,subs_dir{cnt}];
    if ~exist(this_dir,'dir')
        mkdir(this_dir)
    end
end