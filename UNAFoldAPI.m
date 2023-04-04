% Define the RNA sequence
sequence = 'GGGAAAUCC';
[RNAbracket, Energy] = rnafold(sequence)
% % Define the folding parameters
% params = struct('T', 37.0, 'sodium', 1.0, 'magnesium', 0.0);
% 
% % Encode the sequence and parameters as a JSON object
% data = struct('sequence', sequence, 'params', params);
% json = jsonencode(data);
% 
% % Send a POST request to the UNAFold API
% url = 'https://www.unafold.org/cgi-bin/UNAfold.cgi';
% options = weboptions('MediaType', 'application/json', 'RequestMethod', 'post', 'Timeout', 10);
% response = webwrite(url, json, options);
% 
% % Decode the response JSON object
% result = jsondecode(response);
% 
% % Extract the predicted secondary structure and folding energies
% structure = result.structure;
% dG = result.dG;