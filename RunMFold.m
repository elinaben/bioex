
function folding_energy_array = runMFold(sequence, sliding_window_size)

    
    sequence = 'GAGTAGTGGAACTAGGCTATGC';
    command = sprintf('echo "%s" | RNAfold ', sequence);
    [status, result] = system(command);
    disp(result);
end


