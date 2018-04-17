#!/usr/bin/env python
#
# Copyright 2018, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of XXX.
#
# XXX is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# XXX is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with XXX.  If not, see <http://www.gnu.org/licenses/>.
#


import os, sys, random, math
import tensorflow as tf
import numpy as np
from datetime import datetime
from argparse import ArgumentParser, FileType

def weight_variable(shape):
  initial = tf.truncated_normal(shape, stddev=0.1)
  return tf.Variable(initial)

def bias_variable(shape):
  initial = tf.constant(0.1, shape=shape)
  return tf.Variable(initial)


"""
"""
def test_circle():
  num_theta = 1000
  data = []
  for i in range(num_theta):
    theta = 2 * math.pi * i / (num_theta - 1)
    data.append([theta, math.cos(theta), math.sin(theta)])

  feature_len = 1
  random.shuffle(data)
  features = [data[i][:feature_len] for i in range(len(data))]
  labels = [data[i][feature_len:] for i in range(len(data))]

  label_len = len(labels[0])
  assert len(features) > 0 and len(features) == len(labels)

  x = tf.placeholder(tf.float32, [None, feature_len])
  y = tf.placeholder(tf.float32, [None, label_len])

  n = 32
  W_fc1 = weight_variable([feature_len, n])
  b_fc1 = bias_variable([n])
  h_fc1 = tf.nn.relu(tf.matmul(x, W_fc1) + b_fc1)

  for m in range(0):
    W_temp = weight_variable([n, n])
    b_temp = bias_variable([n])
    h_fc1 = tf.nn.relu(tf.matmul(h_fc1, W_temp) + b_temp)

  W_fc2 = weight_variable([n, label_len])
  b_fc2 = bias_variable([label_len])
  y_ = tf.matmul(h_fc1, W_fc2) + b_fc2

  meansq = tf.reduce_mean(tf.square(y - y_))
  l1_regularizer = tf.contrib.layers.l1_regularizer(scale=0.005, scope=None)
  weights = tf.trainable_variables() # all vars of your graph
  regularization_penalty = tf.contrib.layers.apply_regularization(l1_regularizer, weights)
  regularized_loss = meansq + regularization_penalty # this loss needs to be minimized
  optimizer = tf.train.AdamOptimizer(1e-4)

  # train = optimizer.minimize(regularized_loss)
  train = optimizer.minimize(meansq)

  sess = tf.Session()
  init = tf.global_variables_initializer()
  sess.run(init)
  for i in range(30000 + 1):
    sess.run(train, {x: features, y: labels})
    if i % 5000 == 0:
        print ("Iteration:", i)
        # print "\t", sess.run([regularized_loss], {x: features, y: labels})
        print ("\t", sess.run([meansq], {x: features, y: labels}))

    for i in range(8):
      theta = math.pi * i / 4
      prediction = sess.run(y_, {x: [[theta]]})
      print ("At %.2fpi, cos is %.4f (%.4f, true) and sin is %.4f (%.4f)" % \
               (float(i) / 4, prediction[0][0], math.cos(theta), prediction[0][1], math.sin(theta)))

    file = open("sphere_mat_i%d.txt" % iter, 'w')
    file.write("%d\n" % len(weights))
    for weight in weights:
      nparray = weight.eval(sess)
      if weight.shape.ndims == 2:
        dim_str = "%d %d" % (nparray.shape[0], nparray.shape[1])
      else:
        dim_str = "%d" % nparray.shape[0]
        wt_str = [str(wt) for wt in np.ravel(nparray)]
        file.write("%s %s\n" % (dim_str, ' '.join(wt_str)))
    file.close()

  
"""
"""
def test_sphere():
  num_theta, num_theta2 = 500, 500
  data = []
  for i in range(num_theta):
    theta = 2 * math.pi * i / (num_theta - 1)
    for j in range(num_theta2):
      theta2 = 2 * math.pi * j / (num_theta2 - 1)
      data.append([theta, theta2, math.cos(theta) * math.sin(theta2), math.sin(theta) * math.sin(theta2), math.cos(theta2)])
      
  feature_len = 2
  random.shuffle(data)
  features = [data[i][:feature_len] for i in range(len(data))]
  labels = [data[i][feature_len:] for i in range(len(data))]
    
  label_len = len(labels[0])
  assert len(features) > 0 and len(features) == len(labels)
  
  x = tf.placeholder(tf.float32, [None, feature_len])
  y = tf.placeholder(tf.float32, [None, label_len])
    
  n = 64
  W_fc1 = weight_variable([feature_len, n])
  b_fc1 = bias_variable([n])
  h_fc1 = tf.nn.relu(tf.matmul(x, W_fc1) + b_fc1)
  
  for m in range(10):
    W_temp = weight_variable([n, n])
    b_temp = bias_variable([n])
    h_fc1 = tf.nn.relu(tf.matmul(h_fc1, W_temp) + b_temp)
    
  W_fc2 = weight_variable([n, label_len])
  b_fc2 = bias_variable([label_len])
  y_ = tf.matmul(h_fc1, W_fc2) + b_fc2
    
  meansq = tf.reduce_mean(tf.square(y - y_))
  l1_regularizer = tf.contrib.layers.l1_regularizer(scale=0.005, scope=None)
  weights = tf.trainable_variables() # all vars of your graph
  regularization_penalty = tf.contrib.layers.apply_regularization(l1_regularizer, weights)
  regularized_loss = meansq + regularization_penalty # this loss needs to be minimized
  optimizer = tf.train.AdamOptimizer(1e-4)
    
  # train = optimizer.minimize(regularized_loss)
  train = optimizer.minimize(meansq)

  sess = tf.Session()
  init = tf.global_variables_initializer()
  sess.run(init)
  for iter in range(100000 + 1):
    sess.run(train, {x: features, y: labels})
    if iter % 1000 == 0:
      print (str(datetime.now()), "Iteration:", iter)
      # print "\t", sess.run([regularized_loss], {x: features, y: labels})
      print ("\t", sess.run([meansq], {x: features, y: labels}))

      for i in range(4):
        theta = math.pi * i / 4
        for j in range(2):
          theta2 = math.pi * j / 2
          prediction = sess.run(y_, {x: [[theta, theta2]]})
          print ("At %.2fpi and %.2fpi, x is %.3f (%.3f, true), y is %.3f (%.3f), and z is %.3f (%.3f)" % \
                   (float(i) / 4, float(j) / 2, prediction[0][0], math.cos(theta) * math.sin(theta2), prediction[0][1], math.sin(theta) * math.sin(theta2), prediction[0][2], math.cos(theta2)))

      file = open("sphere_mat_i%d.txt" % iter, 'w')
      file.write("%d\n" % len(weights))
      for weight in weights:
        nparray = weight.eval(sess)
        if weight.shape.ndims == 2:
          dim_str = "%d %d" % (nparray.shape[0], nparray.shape[1])
        else:
          dim_str = "%d" % nparray.shape[0]
          wt_str = [str(wt) for wt in np.ravel(nparray)]
          file.write("%s %s\n" % (dim_str, ' '.join(wt_str)))
      file.close()

        
"""
"""
if __name__ == '__main__':
    # test_circle()
    test_sphere()
