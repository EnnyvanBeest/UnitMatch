# Wentao Qiu, 2023-10-07
# qiuwentao1212@gmail.com

import torch
from torch import nn


### Conv module using LayerNorm and GELU
class Conv1DSeq(nn.Module):
    """
    Conv1DSeq layer used to convolve input data.
    """
    def __init__(self, n_channels, n_filters, kernel_sizes, dilation_rates, expected_seq_len):
        """
        Initialize `Conv1DSeq` object.
        :param n_filters: (2[list],) - The dimensionality of the output space.
        :param kernel_sizes: (2[list],) - The length of the 1D convolution window.
        :param dilation_rates: (2[list],) - The dilation rate to use for dilated convolution.
        :param expected_seq_len: int - The expected sequence length after convolution.
        """
        super(Conv1DSeq, self).__init__()

        assert len(n_filters) == len(kernel_sizes) == len(dilation_rates) == 2
        self.n_filters = n_filters
        self.kernel_sizes = kernel_sizes
        self.dilation_rates = dilation_rates

        # Initialize the first component of `Conv1DSeq`.
        self.conv1 = nn.Conv1d(n_channels, n_filters[0], kernel_sizes[0], padding="same",
                               dilation=dilation_rates[0])
        self.ln1 = nn.LayerNorm([n_filters[0],expected_seq_len])
        
        # Initialize the second component of `Conv1DSeq`.
        self.conv2 = nn.Conv1d(n_filters[0], n_filters[1], kernel_sizes[1], padding="same",
                               dilation=dilation_rates[1])
        self.ln2 = nn.LayerNorm([n_filters[1], expected_seq_len])
        self.gelu = nn.GELU()

    def forward(self, inputs):
        """
        :param inputs: (batch_size, seq_len, n_channels) - The input data.
        :return outputs: (batch_size, seq_len, n_output_channels) - The convolved data.
        """
        outputs = self.conv1(inputs.permute(0, 2, 1)) + inputs.permute(0, 2, 1)
        # outputs = nn.functional.gelu(self.ln1(outputs))
        outputs = self.gelu(self.ln1(outputs))
        outputs = self.conv2(outputs) + outputs
        # outputs = nn.functional.gelu(self.ln2(outputs))
        outputs = self.gelu(self.ln2(outputs))
        return outputs
    
class SpatioTemporalCNN_V2(nn.Module):
    def __init__(self,n_channel,n_time,n_output=256):
        super().__init__()
        self.n_channel = n_channel
        self.n_time = n_time
        self.n_output = n_output
        self.ConvBlock1 = Conv1DSeq(self.n_channel, [self.n_channel, self.n_channel],[3,3],[1,2],self.n_time)
        self.ConvBlock2 = Conv1DSeq(self.n_time, [self.n_time, self.n_time],[3,3],[1,2],self.n_channel)
        
        self.n_channel_red = self.n_channel//2
        self.SpatialBlock = nn.Sequential(
            nn.Conv1d(self.n_channel,self.n_channel_red, kernel_size=1),
            nn.LayerNorm([self.n_channel_red, self.n_time]),
            nn.GELU(),
        )

        self.FcBlock = nn.Sequential(
            nn.Linear(self.n_channel_red*self.n_time, 512),
            nn.LayerNorm(512),
            nn.GELU(),
            nn.Linear(512, 512),
            nn.LayerNorm(512),
            nn.GELU(),
            nn.Linear(512, self.n_output),
        )

        self.__init_weight()
    
    def __init_weight(self):
        for m in self.modules():
            if isinstance(m, nn.Conv1d) or isinstance(m, nn.Linear):
                nn.init.kaiming_normal_(m.weight)
                if m.bias is not None:
                    m.bias.data.zero_()
            elif isinstance(m, nn.LayerNorm):
                m.weight.data.fill_(1)
                m.bias.data.zero_()
        
    def forward(self, x):
        """
        :inputs: (batch_size, n_time, n_channels) 
        :return outputs: (batch_size, self.n_output)
        """
        x = self.ConvBlock1(x) # after this, the shape is [bsz, C, T]
        x = self.ConvBlock2(x) # after this, the shape is [bsz, T, C]
        x = x.permute(0,2,1)
        x = self.SpatialBlock(x) # after this, the shape is [bsz, C_red, T]
        x = x.view(x.shape[0],-1)
        x = self.FcBlock(x)
        return x

class Decoder_SpatioTemporalCNN_V2(nn.Module):
    def __init__(self, n_channel,n_time,n_input=256):
        super().__init__()
        self.n_channel = n_channel
        self.n_time = n_time
        self.n_channel_red = self.n_channel//2
        self.FcBlock = nn.Sequential(
            nn.Linear(n_input, 512),
            nn.LayerNorm(512),
            nn.GELU(),
            nn.Linear(512, self.n_channel_red*self.n_time),
        )
        self.DeConvSpatialBlock = nn.Sequential(
            nn.ConvTranspose1d(self.n_channel_red, self.n_channel, kernel_size=1),
            nn.LayerNorm([self.n_channel, self.n_time]),
            nn.GELU(),
        )
        self.DeConvSpatialBlock2 = nn.Sequential(
            nn.ConvTranspose1d(self.n_channel, self.n_channel, kernel_size=3, padding=1),
            nn.LayerNorm([self.n_channel, self.n_time]),
            nn.GELU(),
        )
        self.DeConvTimeLayer = nn.ConvTranspose1d(self.n_time, self.n_time, kernel_size=3, padding=1)
        self.__init_weight()
    
    def __init_weight(self):
        for m in self.modules():
            if isinstance(m, nn.ConvTranspose1d) or isinstance(m, nn.Linear):
                nn.init.kaiming_normal_(m.weight)
                if m.bias is not None:
                    m.bias.data.zero_()
            elif isinstance(m, nn.LayerNorm):
                m.weight.data.fill_(1)
                m.bias.data.zero_()

    def forward(self, x):
        '''
        Inputs: (batch_size, n_input)
        Outputs: (batch_size, n_time, n_channel)
        '''
        x = self.FcBlock(x)
        x = x.view(x.shape[0],self.n_channel_red,self.n_time)
        x = self.DeConvSpatialBlock(x) # after this, the shape is [bsz, C, T]
        x = self.DeConvSpatialBlock2(x)
        x = x.permute(0,2,1)
        x = self.DeConvTimeLayer(x)
        return x

class SpatioTemporalAutoEncoder_V2(nn.Module):
    def __init__(self,n_channel,n_time,n_output=256):
        super().__init__()
        self.n_channel = n_channel
        self.n_time = n_time
        self.n_output = n_output
        self.encoder = SpatioTemporalCNN_V2(self.n_channel,self.n_time,self.n_output)
        self.decoder = Decoder_SpatioTemporalCNN_V2(self.n_channel,self.n_time,self.n_output)

    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded

def count_parameters(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)

if __name__ == '__main__':
    # Create a random input tensor.
    inputs = torch.randn(5, 60, 30)
    
    ### Test SpatioTemporalCNN_V2
    AE_V2_encoder = SpatioTemporalCNN_V2(30, 60, 256)
    AE_V2_decoder = Decoder_SpatioTemporalCNN_V2(30, 60, 256)
    AE_V2 = SpatioTemporalAutoEncoder_V2(30, 60, 256)

    # # Count the number of parameters in the model.
    print("Number of parameters in AE_V2_encoder:", count_parameters(AE_V2_encoder))
    print("Number of parameters in AE_V2_decoder:", count_parameters(AE_V2_decoder))
    print("Number of parameters in SpatioTemporalCNN_V2:", count_parameters(AE_V2))
