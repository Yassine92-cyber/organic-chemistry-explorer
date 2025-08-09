import React from 'react';
import { X, CheckCircle, AlertCircle, Info, Download, FileText } from 'lucide-react';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from '@/components/ui/table';
import { Separator } from '@/components/ui/separator';

export interface ImportError {
  row: number;
  smiles: string;
  error: string;
  suggestion?: string;
}

export interface ImportReport {
  totalRows: number;
  successfulImports: number;
  failedImports: number;
  errors: ImportError[];
  warnings: string[];
  executionTime: number;
  source: string;
}

interface ImportReportModalProps {
  report: ImportReport | null;
  isOpen: boolean;
  onClose: () => void;
}

const ImportReportModal: React.FC<ImportReportModalProps> = ({
  report,
  isOpen,
  onClose
}) => {
  if (!isOpen || !report) return null;

  const successRate = ((report.successfulImports / report.totalRows) * 100).toFixed(1);
  const errorRate = ((report.failedImports / report.totalRows) * 100).toFixed(1);

  const exportReport = () => {
    const reportData = {
      timestamp: new Date().toISOString(),
      ...report,
      summary: {
        successRate: `${successRate}%`,
        errorRate: `${errorRate}%`,
        executionTime: `${report.executionTime}ms`
      }
    };

    const blob = new Blob([JSON.stringify(reportData, null, 2)], {
      type: 'application/json'
    });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `import_report_${new Date().toISOString().split('T')[0]}.json`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  return (
    <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center z-50 p-4">
      <div className="bg-white rounded-lg w-full max-w-4xl max-h-[90vh] flex flex-col">
        <div className="flex items-center justify-between p-6 border-b">
          <div className="flex items-center space-x-3">
            <FileText className="w-6 h-6 text-blue-600" />
            <div>
              <h2 className="text-xl font-bold">Import Report</h2>
              <p className="text-sm text-gray-500">Source: {report.source}</p>
            </div>
          </div>
          <div className="flex items-center space-x-2">
            <Button
              variant="outline"
              size="sm"
              onClick={exportReport}
            >
              <Download className="w-4 h-4 mr-2" />
              Export
            </Button>
            <Button
              variant="ghost"
              size="sm"
              onClick={onClose}
            >
              <X className="w-4 h-4" />
            </Button>
          </div>
        </div>

        <div className="flex-1 overflow-auto p-6 space-y-6">
          {/* Summary Cards */}
          <div className="grid grid-cols-1 md:grid-cols-4 gap-4">
            <Card>
              <CardContent className="p-4">
                <div className="flex items-center space-x-2">
                  <div className="w-3 h-3 bg-blue-500 rounded-full"></div>
                  <span className="text-sm font-medium">Total Rows</span>
                </div>
                <p className="text-2xl font-bold mt-2">{report.totalRows}</p>
              </CardContent>
            </Card>

            <Card>
              <CardContent className="p-4">
                <div className="flex items-center space-x-2">
                  <CheckCircle className="w-4 h-4 text-green-600" />
                  <span className="text-sm font-medium">Successful</span>
                </div>
                <p className="text-2xl font-bold mt-2 text-green-600">
                  {report.successfulImports}
                </p>
                <p className="text-xs text-gray-500">{successRate}%</p>
              </CardContent>
            </Card>

            <Card>
              <CardContent className="p-4">
                <div className="flex items-center space-x-2">
                  <AlertCircle className="w-4 h-4 text-red-600" />
                  <span className="text-sm font-medium">Failed</span>
                </div>
                <p className="text-2xl font-bold mt-2 text-red-600">
                  {report.failedImports}
                </p>
                <p className="text-xs text-gray-500">{errorRate}%</p>
              </CardContent>
            </Card>

            <Card>
              <CardContent className="p-4">
                <div className="flex items-center space-x-2">
                  <Info className="w-4 h-4 text-blue-600" />
                  <span className="text-sm font-medium">Time</span>
                </div>
                <p className="text-2xl font-bold mt-2">
                  {report.executionTime}ms
                </p>
              </CardContent>
            </Card>
          </div>

          {/* Warnings */}
          {report.warnings.length > 0 && (
            <Card>
              <CardHeader>
                <CardTitle className="text-lg flex items-center space-x-2">
                  <AlertCircle className="w-5 h-5 text-yellow-600" />
                  <span>Warnings</span>
                  <Badge variant="secondary">{report.warnings.length}</Badge>
                </CardTitle>
              </CardHeader>
              <CardContent>
                <ul className="space-y-2">
                  {report.warnings.map((warning, index) => (
                    <li key={index} className="flex items-start space-x-2 text-sm">
                      <span className="text-yellow-600 mt-0.5">•</span>
                      <span className="text-gray-700">{warning}</span>
                    </li>
                  ))}
                </ul>
              </CardContent>
            </Card>
          )}

          {/* Error Details */}
          {report.errors.length > 0 && (
            <Card>
              <CardHeader>
                <CardTitle className="text-lg flex items-center space-x-2">
                  <AlertCircle className="w-5 h-5 text-red-600" />
                  <span>Error Details</span>
                  <Badge variant="destructive">{report.errors.length}</Badge>
                </CardTitle>
              </CardHeader>
              <CardContent>
                <div className="max-h-96 overflow-auto">
                  <Table>
                    <TableHeader>
                      <TableRow>
                        <TableHead className="w-16">Row</TableHead>
                        <TableHead className="w-48">SMILES</TableHead>
                        <TableHead>Error</TableHead>
                        <TableHead className="w-48">Suggestion</TableHead>
                      </TableRow>
                    </TableHeader>
                    <TableBody>
                      {report.errors.map((error, index) => (
                        <TableRow key={index}>
                          <TableCell className="font-mono text-sm">
                            {error.row}
                          </TableCell>
                          <TableCell className="font-mono text-sm max-w-48 truncate" title={error.smiles}>
                            {error.smiles}
                          </TableCell>
                          <TableCell className="text-sm text-red-600">
                            {error.error}
                          </TableCell>
                          <TableCell className="text-sm text-gray-600">
                            {error.suggestion || '-'}
                          </TableCell>
                        </TableRow>
                      ))}
                    </TableBody>
                  </Table>
                </div>
              </CardContent>
            </Card>
          )}
        </div>

        <div className="p-6 border-t bg-gray-50">
          <div className="flex justify-between items-center">
            <p className="text-sm text-gray-600">
              Import completed at {new Date().toLocaleTimeString()}
            </p>
            <Button onClick={onClose}>
              Close
            </Button>
          </div>
        </div>
      </div>
    </div>
  );
};

export default ImportReportModal; 