import logging
import os, sys
import argparse
import subprocess
import numpy as np
import tqdm 
import torch
import torch.optim as optim
from torch.utils.data import DataLoader
from torch.utils.tensorboard  import SummaryWriter
from pathlib import Path

from utils import metric
from utils.losses import *
from utils.npdataset import NeuropixelsDataset_cortexlab, TrainExperimentBatchSampler, ValidationExperimentBatchSampler
from utils.mymodel import *

logger = logging.getLogger(__name__)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def validation(epoch, model, projector, val_loader, clip_loss, writer):
    model.eval()
    clip_loss.eval()
    losses = metric.AverageMeter()
    experiment_accuracies = []
    if torch.cuda.is_available():
        model = model.cuda()
        projector = projector.cuda()
        clip_loss = clip_loss.cuda()

    with torch.no_grad():
        progress_bar = tqdm.tqdm(total=len(val_loader), desc='Epoch {:3d}'.format(epoch))
        for estimates, candidates,_,_,_ in val_loader:
            if torch.cuda.is_available():
                estimates = estimates.cuda()
                candidates = candidates.cuda()

            bsz = estimates.shape[0]
            # Forward pass
            enc_estimates = model(estimates) # shape [bsz, channel*time]
            enc_candidates = model(candidates) # shape [bsz, channel*time]
            proj_estimates = projector(enc_estimates)
            proj_candidates = projector(enc_candidates)
            loss_clip = clip_loss(proj_estimates, proj_candidates)
            loss = loss_clip
            losses.update(loss.item(), bsz)

            probs = clip_prob(enc_estimates, enc_candidates)
            predicted_indices = torch.argmax(probs, dim=1)  # Get the index of the max probability for each batch element
            ground_truth_indices = torch.arange(bsz, device=device)  # Diagonal indices as ground truth
            correct_predictions = (predicted_indices == ground_truth_indices).sum().item()  # Count correct predictions
            accuracy = correct_predictions / bsz
            experiment_accuracies.append(accuracy)
            progress_bar.update(1)
        progress_bar.close()
    
    print('Epoch: %d'%(epoch), 'Validation Loss: %.9f'%(losses.avg), 'Validation Accuracy: %.9f'%(np.mean(experiment_accuracies)))
    # print('clip temp tau', clip_loss.temp_tau)
    writer.add_scalar('Validation/Loss', losses.avg, epoch)
    writer.add_scalar('Validation/Accuracy', np.mean(experiment_accuracies), epoch)
    return

def train(epoch, model, projector, optimizer, train_loader, clip_loss, writer):
    model.train()
    clip_loss.train()
    losses = metric.AverageMeter()
    iteration = len(train_loader) * epoch
    if torch.cuda.is_available():
        model = model.cuda()
        projector = projector.cuda()
        clip_loss = clip_loss.cuda()

    progress_bar = tqdm.tqdm(total=len(train_loader), desc='Epoch {:3d}'.format(epoch))
    for estimates, candidates,_,_,_ in train_loader:
        bsz = estimates.shape[0]
        if torch.cuda.is_available():
            estimates = estimates.cuda()
            candidates = candidates.cuda()
        optimizer.zero_grad()
        enc_estimates = model(estimates) # shape [bsz, channel*time]
        enc_candidates = model(candidates) # shape [bsz, channel*time]
        proj_estimates = projector(enc_estimates)
        proj_candidates = projector(enc_candidates)
        loss_clip = clip_loss(proj_estimates, proj_candidates)
        loss = loss_clip
        losses.update(loss.item(), bsz)
        loss.backward()
        optimizer.step()
        progress_bar.update(1)
        iteration += 1
        if iteration % 50 == 0:
            writer.add_scalar('Train/Loss', losses.avg, iteration)
    
    progress_bar.close()
    print(' Epoch: %d'%(epoch), 'Loss: %.9f'%(losses.avg))
    return 

def run(args):
    save_folder = os.path.join('ModelExp','experiments',  args.exp_name)
    ckpt_folder = os.path.join(save_folder,  'ckpt')
    log_folder = os.path.join(save_folder,  'log')
    os.makedirs(ckpt_folder, exist_ok=True)
    os.makedirs(log_folder, exist_ok=True)
    writer = SummaryWriter(log_dir=log_folder)

    train_data_root = args.train_root
    np_dataset = NeuropixelsDataset_cortexlab(data_dir=train_data_root, batch_size=args.batchsize, mode='train')
    train_sampler = TrainExperimentBatchSampler(np_dataset, args.batchsize, shuffle=True)
    train_loader = DataLoader(np_dataset, batch_sampler=train_sampler)
    val_np_dataset = NeuropixelsDataset_cortexlab(data_dir=train_data_root, batch_size=args.batchsize, mode='train')
    val_sampler = ValidationExperimentBatchSampler(val_np_dataset, shuffle = True)
    val_loader = DataLoader(val_np_dataset, batch_sampler=val_sampler)

    print('train dataset length: %d'%(len(np_dataset)))

    print(f"To open Tensorboard, run this: tensorboard --logdir {os.path.join(os.getcwd(), log_folder)}")

    model = SpatioTemporalCNN_V2(n_channel=30,n_time=60,n_output=256).to(device)
    
    model = model.double()
    finetune_folder = os.path.join('ModelExp','AE_experiments',  args.finetune)
    ckpt_finetune_folder = os.path.join(finetune_folder,  'ckpt')
    ckpt_lst = os.listdir(ckpt_finetune_folder)
    ckpt_lst.sort(key=lambda x: int(x.split('_')[-1]))
    read_path = os.path.join(ckpt_finetune_folder, ckpt_lst[-1])
    print('load checkpoint from %s'%(read_path))
    checkpoint = torch.load(read_path)
    model.load_state_dict(checkpoint['encoder'])
    for name, param in model.named_parameters():
        if 'FcBlock' not in name:
            param.requires_grad = False

    projector = Projector(input_dim=256, output_dim=128, hidden_dim=128, n_hidden_layers=1, dropout=0.1).to(device)
    projector = projector.double()

    # clip_loss = ClipLoss1D().to(device)
    clip_loss = CustomClipLoss().to(device)

    encoder_fc_params = [param for name, param in model.named_parameters() if 'FcBlock' in name and param.requires_grad]
    projector_params = list(projector.parameters())  # Assuming projector is defined elsewhere
    clip_loss_params = list(clip_loss.parameters())  # Assuming clip_loss is defined elsewhere

    # Combine parameters from different parts with their respective learning rates
    optimizer_params = [
        {'params': encoder_fc_params, 'lr': args.lr_enc},  # Smaller learning rate for FcBlock
        {'params': projector_params + clip_loss_params, 'lr': args.lr_proj}  # Larger learning rate for projector and clip_loss
    ]

    optimizer = optim.Adam(optimizer_params)

    if args.cont:
        # load latest checkpoint
        ckpt_lst = os.listdir(ckpt_folder)
        ckpt_lst.sort(key=lambda x: int(x.split('_')[-1]))
        read_path = os.path.join(ckpt_folder, ckpt_lst[-1])
        print('load checkpoint from %s'%(read_path))
        checkpoint = torch.load(read_path)
        model.load_state_dict(checkpoint['model'])
        optimizer.load_state_dict(checkpoint['optimizer'])
        clip_loss.load_state_dict(checkpoint['clip_loss'])
        start_epoch = checkpoint['epoch'] + 1
    else:
        start_epoch = 0

    if args.total_epoch == 0:
        # don't train, just save the untrained checkpoint
        state = {
            'model': model.state_dict(),
            'optimizer': optimizer.state_dict(),
            'clip_loss': clip_loss.state_dict(),
            'epoch': 0,
        }
        save_file = os.path.join(ckpt_folder, 'ckpt_epoch_0')
        torch.save(state, save_file)

    for epoch in range(start_epoch, args.total_epoch):
        train(epoch, model, projector, optimizer, train_loader, clip_loss, writer)
        if epoch % args.save_freq == 0:
            state = {
                'model': model.state_dict(),
                'optimizer': optimizer.state_dict(),
                'clip_loss': clip_loss.state_dict(),
                'epoch': epoch,
            }
            save_file = os.path.join(ckpt_folder, 'ckpt_epoch_%s'%(str(epoch)))
            torch.save(state, save_file)

        # validate and test
        validation(epoch, model, projector, val_loader, clip_loss, writer)
    # test(epoch, model, test_loader, writer)
    
    return

def run_finetune(exp_name, dataset, lr_enc=2*1e-5, lr_proj=1.1*1e-4, save_freq=1, total_epoch=50, cont=False, batchsize=40):

    from argparse import Namespace
    args = Namespace(
        exp_name=exp_name,
        lr_enc=lr_enc,
        lr_proj=lr_proj,
        save_freq=save_freq,
        total_epoch=total_epoch,
        cont=cont,
        batchsize=batchsize,
        dataset=dataset
    )

    save_path = Path(__file__).parent.parent
    save_folder = os.path.join(save_path, 'ModelExp', 'experiments', args.exp_name)
    if os.path.exists(save_folder) and not args.cont:
        # Check if the folder already exists and if we are not continuing training
        raise ValueError("Running this experiment will overwrite a previous one. Please choose a different name.")
    
    ckpt_folder = os.path.join(save_folder,  'ckpt')
    log_folder = os.path.join(save_folder,  'log')
    os.makedirs(ckpt_folder, exist_ok=True)
    os.makedirs(log_folder, exist_ok=True)
    writer = SummaryWriter(log_dir=log_folder)

    np_dataset = dataset
    train_sampler = TrainExperimentBatchSampler(np_dataset, args.batchsize, shuffle=True)
    train_loader = DataLoader(np_dataset, batch_sampler=train_sampler)

    val_sampler = ValidationExperimentBatchSampler(dataset, shuffle=True)
    val_loader = DataLoader(dataset, batch_sampler=val_sampler)

    print('train dataset length: %d'%(len(np_dataset)))

    import psutil
    # Kill any existing process on port 6006
    for proc in psutil.process_iter(['pid', 'name']):
        try:
            if proc.name() == 'tensorboard.exe':
                proc.kill()
        except:
            pass

    tensorboard_process = subprocess.Popen(['tensorboard', '--logdir', log_folder, '--port', '6006'],
                                           stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE)
    print(f"To view the tensorboard, click this link: http://localhost:6006")

    model = SpatioTemporalCNN_V2(n_channel=30,n_time=60,n_output=256).to(device)
    
    model = model.double()
    finetune_folder = os.path.join(save_path, 'ModelExp','AE_experiments',  args.exp_name)
    ckpt_finetune_folder = os.path.join(finetune_folder,  'ckpt')
    ckpt_lst = os.listdir(ckpt_finetune_folder)
    ckpt_lst.sort(key=lambda x: int(x.split('_')[-1]))
    read_path = os.path.join(ckpt_finetune_folder, ckpt_lst[-1])
    print('load checkpoint from %s'%(read_path))
    checkpoint = torch.load(read_path)
    model.load_state_dict(checkpoint['encoder'])
    for name, param in model.named_parameters():
        if 'FcBlock' not in name:
            param.requires_grad = False

    projector = Projector(input_dim=256, output_dim=128, hidden_dim=128, n_hidden_layers=1, dropout=0.1).to(device)
    projector = projector.double()

    clip_loss = CustomClipLoss().to(device)

    encoder_fc_params = [param for name, param in model.named_parameters() if 'FcBlock' in name and param.requires_grad]
    projector_params = list(projector.parameters())
    clip_loss_params = list(clip_loss.parameters())

    # Combine parameters from different parts with their respective learning rates
    optimizer_params = [
        {'params': encoder_fc_params, 'lr': args.lr_enc},  # Smaller learning rate for FcBlock
        {'params': projector_params + clip_loss_params, 'lr': args.lr_proj}  # Larger learning rate for projector and clip_loss
    ]

    optimizer = optim.Adam(optimizer_params)

    if args.cont:
        # load latest checkpoint
        ckpt_lst = os.listdir(ckpt_folder)
        ckpt_lst.sort(key=lambda x: int(x.split('_')[-1]))
        read_path = os.path.join(ckpt_folder, ckpt_lst[-1])
        print('load checkpoint from %s'%(read_path))
        checkpoint = torch.load(read_path)
        model.load_state_dict(checkpoint['model'])
        optimizer.load_state_dict(checkpoint['optimizer'])
        clip_loss.load_state_dict(checkpoint['clip_loss'])
        start_epoch = checkpoint['epoch'] + 1
    else:
        start_epoch = 0

    if args.total_epoch == 0:
        # don't train, just save the untrained checkpoint
        state = {
            'model': model.state_dict(),
            'optimizer': optimizer.state_dict(),
            'clip_loss': clip_loss.state_dict(),
            'epoch': 0,
        }
        save_file = os.path.join(ckpt_folder, 'ckpt_epoch_0')
        torch.save(state, save_file)

    for epoch in range(start_epoch, args.total_epoch):
        train(epoch, model, projector, optimizer, train_loader, clip_loss, writer)
        if epoch % args.save_freq == 0:
            state = {
                'model': model.state_dict(),
                'optimizer': optimizer.state_dict(),
                'clip_loss': clip_loss.state_dict(),
                'epoch': epoch,
            }
            save_file = os.path.join(ckpt_folder, 'ckpt_epoch_%s'%(str(epoch)))
            torch.save(state, save_file)

        # validate
        validation(epoch, model, projector, val_loader, clip_loss, writer)
    
    return

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('--exp_name', '-e', type=str, required=True, help="The checkpoints and logs will be save in ./checkpoint/$EXP_NAME")
    arg_parser.add_argument('--finetune', '-f', type=str, required=True, help="Load the AE encoder from the path ./checkpoint/$finetune")
    arg_parser.add_argument('--lr_enc', '-le', type=float, default=2*1e-5, help="Learning rate for encoder")
    arg_parser.add_argument('--lr_proj', '-lp', type=float, default=1.1*1e-4, help="Learning rate for projector")
    arg_parser.add_argument('--save_freq', '-s', type=int, default=1, help="frequency of saving model")
    arg_parser.add_argument('--total_epoch', '-t', type=int, default=50, help="total epoch number for training")
    arg_parser.add_argument('--cont', '-c', action='store_true', help="whether to load saved checkpoints from $EXP_NAME and continue training")
    arg_parser.add_argument('--batchsize', '-b', type=int, default=40, help="batch size")
    arg_parser.add_argument('--train_root', type=str, default=r"\path\to\your\data", help="root directory of training data")
    args = arg_parser.parse_args()

    run(args)