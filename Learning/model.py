import torch
import torch.nn as nn

class DeepDrug3D(nn.Module):
    def __init__(self, in_channel):
        super(DeepDrug3D, self).__init__()

        self.conv1 = nn.Conv3d(in_channel, 64, 5)
        self.conv2 = nn.Conv3d(64, 64, 3)

        self.pool = nn.MaxPool3d((2,2,2), stride=None)
        self.fc1 = nn.Linear(64*13*13*13, 128)
        self.fc2 = nn.Linear(128,3)

    def reset_parameters(self):
        self.conv1.reset_parameters()
        self.conv2.reset_parameters()
        self.fc1.reset_parameters()
        self.fc2.reset_parameters()

    def forward(self, x):
        x = self.conv1(x)
        x = nn.LeakyReLU(negative_slope=0.1)(x)
        x = nn.Dropout(p=0.2)(x)

        x = self.conv2(x)
        x = nn.LeakyReLU(negative_slope=0.1)(x)

        x = self.pool(x)
        x = nn.Dropout(p=0.4)(x)

        x = torch.flatten(x, start_dim=1)
        x = self.fc1(x)
        x = nn.LeakyReLU(negative_slope=0.1)(x)
        x = nn.Dropout(p=0.4)(x)

        x = self.fc2(x)
        return x