﻿using System;
using System.IO;
using NLog;

namespace BenchTool
{
    public class TempFolder : IDisposable
    {
        private static Logger Logger { get; } = LogManager.GetCurrentClassLogger();

        public string FullPath { get; private set; }

        public string RootDirName { get; set; }

        public TempFolder(string identifier = "")
        {
            RootDirName = $"{identifier}_{Path.GetFileNameWithoutExtension(Path.GetRandomFileName())}";
            FullPath = Path.Combine(Path.GetTempPath(), RootDirName);
        }

        public void CreateIfNotExist()
        {
            if (!FullPath.IsEmptyOrWhiteSpace())
            {
                FullPath.CreateDirectoryIfNotExist();
            }
        }

        public void Dispose()
        {
            if (Directory.Exists(FullPath))
            {
                try
                {
                    Directory.Delete(FullPath, recursive: true);
                }
                catch (IOException e)
                {
                    Logger.Warn($"{e.Message} {FullPath}");
                }
            }

            FullPath = null;
        }
    }

    public class TempFile : IDisposable
    {
        private static Logger Logger { get; } = LogManager.GetCurrentClassLogger();

        public string FullPath { get; private set; }
        public TempFile()
        {
            string fileName = Path.GetFileNameWithoutExtension(Path.GetRandomFileName());
            //FullPath = Path.Combine(Environment.CurrentDirectory, ".tmp", fileName);
            FullPath = Path.Combine(Path.GetTempPath(), fileName);
        }

        public void Dispose()
        {
            if (File.Exists(FullPath))
            {
                try
                {
                    File.Delete(FullPath);
                }
                catch (IOException e)
                {
                    Logger.Warn($"{e.Message} {FullPath}");
                }
            }

            FullPath = null;
        }
    }
}
